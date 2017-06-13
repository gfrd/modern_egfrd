#ifndef EGFRD_SIMULATOR_HPP
#define EGFRD_SIMULATOR_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "Logger.hpp"
#include "ShellID.hpp"
#include "SerialIDGenerator.hpp"
#include "EventScheduler.hpp"
#include "DomainID.hpp"
#include "World.hpp"
#include "Domain.hpp"
#include "Shell.hpp"
#include "Single.hpp"
#include "ParticleSimulator.hpp"
#include "ShellCreateUtils.hpp"
#include "Pair.hpp"
#include "Multi.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class EGFRDSimulator : public ParticleSimulator
{
public:
   using base_type = ParticleSimulator;
   using boundary_type = World::boundary_type;

   using domain_map = std::unordered_map<DomainID, std::unique_ptr<Domain>>;
   using shell_matrix_type = MatrixSpace<Shell, ShellID, World::particle_matrix_type::SizeX, World::particle_matrix_type::SizeY, World::particle_matrix_type::SizeZ>;
   using particle_id_pair = Particle::particle_id_pair;
   using vec_pid_pair = std::vector<particle_id_pair>;
   using position_structid_pair = Single::position_structid_pair;

   using shell_id_pair = Shell::shell_id_pair;
   using shell_id_pair_and_distance = std::pair<shell_id_pair, double>;
   using shell_id_pair_and_distance_list = std::vector<shell_id_pair_and_distance>;

   // --------------------------------------------------------------------------------------------------------------------------------

   explicit EGFRDSimulator(World& world, ReactionRuleCollection& reaction_rules, RandomNumberGenerator& rng) noexcept
      : base_type(world, reaction_rules, rng), sidgen_(), didgen_(), scheduler_(), shellmat_() { }

   // --------------------------------------------------------------------------------------------------------------------------------

   gi::iteratorRange<domain_map> get_domains() const { return gi::iteratorRange<domain_map>(domains_); }

   // --------------------------------------------------------------------------------------------------------------------------------

   GFRD_EXPORT void dump(std::string filename, bool append = false) const;

   // --------------------------------------------------------------------------------------------------------------------------------

   void add_extrnal_event(double time, CustomAction* action)
   {
      THROW_UNLESS_MSG(illegal_argument, action != nullptr, "No method provided");
      THROW_UNLESS_MSG(illegal_argument, time >= time_, "Cannot insert event in the past");
      scheduler_.add(Event(time, action));
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool step() { return step_intern(); }

   void check() const { check_intern(); }

   void burst_all()
   {
      // create list of domains to burst
      std::vector<DomainID> domains;
      for (const auto&d : domains_)
         domains.emplace_back(d.first);

      while (domains.size() > 0)
      {
         // loop until list is empty, take last
         std::vector<DomainID> ignore;
         DomainID did = domains.back();
         domains.pop_back();

         // if domain is not a Single INIT domain
         const auto& domain = domains_[did];
         auto single = dynamic_cast<Single*>(domain.get());
         if (single != nullptr && single->shell().code() == Shell::Code::INIT) continue;

         // burst it
         burst_domain(did, ignore);

         // remove included bursteds from list
         for (auto id : ignore)
         {
            auto i = std::find(domains.begin(), domains.end(), id);
            if (i != domains.end())
            {
               std::swap(*i, domains.back());
               domains.pop_back();
            }
         }
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void reset()
   {
      scheduler_.clear();
      shellmat_.clear();
      domains_.clear();
      dt_ = time_ = 0;
      num_steps_ = 0;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:

   std::unique_ptr<Single> create_default_single(const DomainID did, const ShellID sid, particle_id_pair pip, std::shared_ptr<Structure> structure, const ReactionRuleCollection::reaction_rule_set& rr) const
   {
      // switch on structure type
      auto *world = dynamic_cast<CuboidalRegion*>(structure.get());
      if (world != nullptr)
      {
         const auto sid_pair = std::make_pair<const ShellID, Shell>(std::move(sid), Shell(did, Sphere(pip.second.position(), pip.second.radius()), Shell::Code::INIT));
         auto domain = std::make_unique<SingleSpherical>(did, pip, sid_pair, rr);
         return std::unique_ptr<Single>(std::move(domain));
      }

      auto *plane = dynamic_cast<PlanarSurface*>(structure.get());
      if (plane != nullptr)
      {
         // TODO scaling parameters , see shells.py:2009
         //auto sid_pair = std::make_pair<const ShellID&, Shell>(sid, Shell(did, Cylinder(pip.second.position(), pip.second.radius(), plane->shape().unit_z(), pip.second.radius()), true));

         // We make a spherical init shell here, construct the cylinder at the next event update
         auto sid_pair = std::make_pair<const ShellID, Shell>(std::move(sid), Shell(did, Sphere(pip.second.position(), pip.second.radius()), Shell::Code::INIT));
         auto domain = std::make_unique<SingleCylindrical>(did, pip, sid_pair, rr);
         return std::unique_ptr<Single>(std::move(domain));
      }

      /*
      elif isinstance(structure, CylindricalSurface):
         testSingle = CylindricalSurfaceSingletestShell(pid_particle_pair, structure, geometrycontainer, domains)
         return CylindricalSurfaceSingle (domain_id, shell_id, testSingle, reaction_rules)
      elif isinstance(structure, DiskSurface):
         try:
            testSingle = CylindricalSurfacePlanarSurfaceInterfaceSingletestShell(pid_particle_pair, structure, geometrycontainer, domains)
         except testShellError as e:
            if __debug__:
               log.warn('Could not make CylindricalSurfacePlanarSurfaceInterfaceSingletestShell, %s' % str(e))
            # making the default testShell should never fail
            testSingle = DiskSurfaceSingletestShell(pid_particle_pair, structure, geometrycontainer, domains)
         return DiskSurfaceSingle (domain_id, shell_id, testSingle, reaction_rules)
      */

      THROW_EXCEPTION(not_implemented, "Structure types not available.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   std::unique_ptr<Pair> create_default_pair(const DomainID did, const ShellID sid, particle_id_pair pip1, std::shared_ptr<Structure> structure1, particle_id_pair pip2, std::shared_ptr<Structure> structure2, const ReactionRuleCollection::reaction_rule_set& rr) const
   {
      if (structure1.get() == structure2.get())
      {
         auto *world = dynamic_cast<CuboidalRegion*>(structure1.get());
         if (world != nullptr)
         {
            auto sid_pair = std::make_pair<const ShellID, Shell>(std::move(sid), Shell(did, Sphere(), Shell::Code::INIT));
            auto domain = std::make_unique<PairSpherical>(did, pip1, pip2, sid_pair, rr);
            return std::unique_ptr<Pair>(std::move(domain));
         }


         //if isinstance(single1.structure, CuboidalRegion) :
         //   return SphericalPairtestShell(single1, single2, geometrycontainer, domains)
         //   elif isinstance(single1.structure, PlanarSurface) :
         //   return PlanarSurfacePairtestShell(single1, single2, geometrycontainer, domains)
         //   elif isinstance(single1.structure, CylindricalSurface) :
         //   return CylindricalSurfacePairtestShell(single1, single2, geometrycontainer, domains)


      }
      else
      {
         // All sorts of Pair Interactions (on different structures)

         // PlanarSurface & PlanarSurface => PlanarSurfaceTransitionPairtestShell
         // PlanarSurface & CuboidalRegion => MixedPair2D3DtestShell
         // CuboidalRegion & PlanarSurface => MixedPair2D3DtestShell
         // PlanarSurface & DiskSurface => MixedPair2DStatictestShell
         // DiskSurface & PlanarSurface  => MixedPair2DStatictestShell
         // CylindricalSurface & DiskSurface => MixedPair1DStatictestShell
         // DiskSurface & CylindricalSurface => MixedPair1DStatictestShell
      }

      THROW_EXCEPTION(not_implemented, "Structure types not available.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   DomainID create_single(particle_id_pair pip)
   {
      // 1. generate identifiers for the domain and shell.
      DomainID did = didgen_();
      ShellID sid = sidgen_();

      const auto& rr = reaction_rules_.query_reaction_rules(pip.second.sid());
      auto species = world_.get_species(pip.second.sid());
      auto structure = world_.get_structure(pip.second.structure_id());

      // 2. Create and register the single domain.The type of the single that will be created depends on the 
      // structure(region or surface) this particle is in / on.Either SphericalSingle, PlanarSurfaceSingle, or CylindricalSurfaceSingle.
      auto single = create_default_single(did, sid, pip, structure, rr);

      // 3. update time and shell container
      single->set_last_time(time_);
      shellmat_.update(single->shell_pair());
      domains_[did] = std::unique_ptr<Domain>(std::move(single));
      return did;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void add_domain_event(DomainID did)
   {
      auto& domain = domains_[did];
      double event_time = time_ + domain->dt();
      EventID eid = scheduler_.add(Event(event_time, did));
      domain->set_event_id(eid);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   const Domain& get_domain(DomainID did) { auto& domain = domains_[did]; return *domain; }

   // --------------------------------------------------------------------------------------------------------------------------------

   void pre_run_check() const
   {
      // check overlaping (or touching) particles at start (user input error)
      for (const auto& p : world_.get_particles())
      {
         auto tpo = world_.test_particle_overlap(p.second.shape(), p.first);
         if (tpo)
         {
            std::string msg;
            msg = make_string() << "Overlapping (or touching) particles in initial configuration: " << static_cast<idtype>(p.first);
            auto ovl = world_.check_particle_overlap(p.second.shape(), p.first);
            for (auto& p2 : ovl) msg = make_string() << msg << " and " << static_cast<idtype>(p2.first.first);
            THROW_EXCEPTION(illegal_argument, msg);
         }
      }

      // particle radius times MULTI_SHELL_FACTOR may not exceed half the matrix cell size (test shell overlap condition)
      double max_particle_radius = shellmat_.cell_size() / (2.0 *GfrdCfg.MULTI_SHELL_FACTOR);
      for (const auto& s : world_.get_species())
      {
         THROW_UNLESS_MSG(illegal_argument, s.second.radius() <= max_particle_radius, "Species '" << static_cast<idtype>(s.first) << ":" << s.second.name() << "' radius (" << s.second.radius() << ") exceeds maximum (" << max_particle_radius << ") for this matrix space.");
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------
   
   void initialize()
   {
      shellmat_.initialize(world_.cell_size());

      pre_run_check();

      set_repulsive();

      // create init-domains for all particles
      for (auto& pip : world_.get_particles())
      {
         DomainID did = create_single(pip);
         add_domain_event(did);
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool step_intern()
   {
      try
      {
         if (num_steps_ == 0) initialize();

         if (scheduler_.empty()) return false;  // nothing to do, simulation done

         auto eid_pair = scheduler_.pop();
         time_ = eid_pair.second.time();
         switch (eid_pair.second.action())
         {
         case Event::actionType::DomainUpdate:
         {
            auto did = eid_pair.second.id();
            dispatch_domain_event(did);
         }  break;

         case Event::actionType::CustomAction:
         {
            auto custom_action = eid_pair.second.custom_action();
            if (custom_action != nullptr)    // Persistence restored simulations may encounter cleared CustomActions.
            {
               custom_action->do_action(time_);
               add_extrnal_event(time_ + custom_action->interval(), custom_action);
            }
         }  break;

         case Event::actionType::Unspecified:
         default:
            THROW_EXCEPTION(illegal_state, "Unknown Event")
         }

         ++num_steps_;

#if defined(_DEBUG)
         //check_intern();
#endif
         return true;
      }
      catch (gfrd_exception ex)
      {
         Log("GFRD").fatal() << ex.what();
         dump(make_string() << "sim_state_" << std::setfill('0') << std::setw(10) << num_steps_ << "_exception" << ".log");
         throw;
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void dispatch_domain_event(DomainID did)
   {
      auto& domain = domains_[did];
      std::vector<DomainID> ignore;

      // do dispatcher thing here (type_info)
      {
         // SINGLE
         auto pse = dynamic_cast<Single*>(domain.get());
         if (pse != nullptr)
         {
            process_single_event(*pse, ignore);
            return;
         }

         // PAIR
         auto ppe = dynamic_cast<Pair*>(domain.get());
         if (ppe != nullptr)
         {
            process_pair_event(*ppe, ignore);
            return;
         }

         // MULTI
         auto pme = dynamic_cast<Multi*>(domain.get());
         if (pme != nullptr)
         {
            process_multi_event(*pme, ignore);
            return;
         }

         THROW_EXCEPTION(illegal_state, "Unknown event");
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   std::vector<particle_id_pair> fire_single_reaction(const particle_id_pair& reactant, const ReactionRule& rrule, const position_structid_pair& pair, std::vector<DomainID>& ignore)
   {
      /*
      # This takes care of the identity change when a single particle decays into one or a number of other particles
      # It performs the reactions:
      # A(any structure) -> 0
      # A(any structure) -> B(same structure or 3D)
      # A(any structure) -> B(same structure) + C(same structure or 3D)
      */

      switch (rrule.get_products().size())
      {
      case 0: // decay to zero
         world_.remove_particle(reactant.first);
         if (rrec_) rrec_->StoreDecayReaction(time_, rrule.id(), reactant.first);
         return std::vector<particle_id_pair>();

      case 1:
      {
         // reactant -> single product
         const SpeciesType& reactant_species = world_.get_species(reactant.second.sid());
         const SpeciesType& product_species = world_.get_species(rrule.get_products()[0]);

         std::vector<Vector3> product_pos_list;
         StructureID product_structure_id;
         if (reactant_species.structure_type_id() != product_species.structure_type_id())
         {
            THROW_EXCEPTION(not_implemented, "Structure change not available.");
         }
         else
         {
            // no change in position
            product_pos_list.emplace_back(reactant.second.position());
            product_structure_id = world_.get_structure(reactant.second.structure_id())->id();
         }

         // make room, when particle relocates or radius is increased
         if (product_species.radius() > reactant_species.radius() || product_pos_list[0] != reactant.second.position())
         {
            Vector3 pos = world_.apply_boundary(product_pos_list[0]);
            double radius = world_.distance(pos, reactant.second.position()) + product_species.radius();
            burst_volume(Sphere(pos, radius * GfrdCfg.SINGLE_SHELL_FACTOR * GfrdCfg.SAFETY), ignore);
         }

         for (auto pos : product_pos_list)
         {
            pos = world_.apply_boundary(pos);
            auto ol = world_.test_particle_overlap(Sphere(pos, product_species.radius()), reactant.first);
            if (!ol)     // no overlap
            {
               world_.remove_particle(reactant.first);
               auto newparticle = world_.add_particle(product_species.id(), product_structure_id, pos);
               if (rrec_) rrec_->StoreUniMolecularReaction(time_, rrule.id(), reactant.first, newparticle.first);
               return std::vector<particle_id_pair>({ newparticle });
            }
         }

         // Process (the lack of) change
         Log("GFRD").warn() << "single reaction: placing product failed, moved reactant instead.";
         auto moved_particle = move_particle(reactant, pair);
         return std::vector<particle_id_pair>({ moved_particle });
      }

      case 2:
      {
         // reactant -> two products
         const SpeciesType& reactant_species = world_.get_species(reactant.second.sid());
         const SpeciesType& product1_species = world_.get_species(rrule.get_products()[0]);
         const SpeciesType& product2_species = world_.get_species(rrule.get_products()[1]);
         double r01 = product1_species.radius() + product2_species.radius();

         if (reactant_species.structure_type_id() != product1_species.structure_type_id() || product1_species.structure_type_id() != product2_species.structure_type_id())
         {
            THROW_EXCEPTION(not_implemented, "Structure change not available.");
         }

         // when fire_single-reaction is called from within a Pair domain, the pair-partner particel may be in the way of two products, so try 4 random placements
         for (int retry = 4; retry > 0; --retry)
         {
            // calculate a random vector in the structure with unit length
            auto reactant_structure = world_.get_structure(reactant.second.structure_id());
            Vector3 iv = GfrdCfg.MINIMAL_SEPARATION_FACTOR * reactant_structure->random_vector(r01, rng_);

            double Dtot = product1_species.D() + product2_species.D();
            Vector3 pos1 = reactant.second.position() - iv * (Dtot != 0.0 ? product1_species.D() / Dtot : 0.5);
            Vector3 pos2 = reactant.second.position() + iv * (Dtot != 0.0 ? product2_species.D() / Dtot : 0.5);

            StructureID product1_structure_id = reactant_structure->id();     // no change in structure
            StructureID product2_structure_id = reactant_structure->id();

            pos1 = world_.apply_boundary(pos1);
            pos2 = world_.apply_boundary(pos2);

            auto ol1 = world_.test_particle_overlap(Sphere(pos1, product1_species.radius()), reactant.first);
            auto ol2 = world_.test_particle_overlap(Sphere(pos2, product2_species.radius()), reactant.first);
            if (!ol1 && !ol2)     // no overlap
            {
               double radius = std::max(world_.distance(reactant.second.position(), pos1) + product1_species.radius() * GfrdCfg.SINGLE_SHELL_FACTOR * GfrdCfg.SAFETY,
                  world_.distance(reactant.second.position(), pos2) + product2_species.radius() * GfrdCfg.SINGLE_SHELL_FACTOR * GfrdCfg.SAFETY);
               burst_volume(Sphere(reactant.second.position(), radius), ignore);    // includes both particles

               world_.remove_particle(reactant.first);
               auto newparticle1 = world_.add_particle(product1_species.id(), product1_structure_id, pos1);
               auto newparticle2 = world_.add_particle(product2_species.id(), product2_structure_id, pos2);
               if (rrec_) rrec_->StoreUnbindingReaction(time_, rrule.id(), reactant.first, newparticle1.first, newparticle2.first);
               return std::vector<particle_id_pair>({ newparticle1, newparticle2 });
            }
            //Log("GFRD").info("single reaction %d: placing product(s) failed, retry %d, iv: %g,%g,%g l=%g", reactant.first, retry, iv.X(), iv.Y(), iv.Z(), iv.length());
         }

         // Process (the lack of) change
         Log("GFRD").warn() << "single reaction " << reactant.first << " placing product(s) failed, moved reactant instead.";
         auto moved_particle = fire_move(reactant, pair);
         return std::vector<particle_id_pair>({ moved_particle });
      }

      default:
         THROW_EXCEPTION(not_implemented, "num products >= 3 not supported.");
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   std::vector<particle_id_pair> fire_pair_reaction(const particle_id_pair& reactant1, const position_structid_pair& pos_pair1,
      const particle_id_pair& reactant2, const position_structid_pair& pos_pair2, const ReactionRule& rrule, Vector3 new_com, std::vector<DomainID>& ignore) const
   {
      UNUSED(ignore);

      // This takes care of the identity change when two particles react with each other
      // It performs the reactions :
      // A(any structure) + B(same structure) -> 0
      // A(any structure) + B(same structure or 3D)->C(same structure)

      switch (rrule.get_products().size())
      {
      case 0: // decay to zero
         world_.remove_particle(reactant1.first);
         world_.remove_particle(reactant2.first);
         if (rrec_) rrec_->StoreAnnihilationReaction(time_, rrule.id(), reactant1.first, reactant2.first);
         return std::vector<particle_id_pair>();

      case 1:
      {
         // reactants -> single product
         //const SpeciesType& reactant1_species = world_.get_species(reactant1.second.sid());
         //const SpeciesType& reactant2_species = world_.get_species(reactant2.second.sid());
         const SpeciesType& product_species = world_.get_species(rrule.get_products()[0]);

         std::vector<Vector3> product_pos_list;
         StructureID product_structure_id;

         // select the structure_id for the product particle to live on.
         // TODO doesn't work for all needed cases
         auto product_structure_ids = world_.get_structure_ids(product_species.structure_type_id());
         if (product_structure_ids.find(pos_pair1.second) != end(product_structure_ids))
            product_structure_id = pos_pair1.second;
         else if (product_structure_ids.find(pos_pair2.second) != end(product_structure_ids))
            product_structure_id = pos_pair2.second;
         else
            THROW_EXCEPTION(unsupported, "reaction product has structureType that is not equal to either one of the reactants.");

         // 3. check that there is space for the products ignoring the reactants
        // Note that we do not check for interfering surfaces(we assume this is no problem)
         auto ol = world_.test_particle_overlap(Sphere(new_com, product_species.radius()), std::vector<ParticleID>{ reactant1.first, reactant2.first});
         if (!ol)     // no overlap
         {
            world_.remove_particle(reactant1.first);
            world_.remove_particle(reactant2.first);
            auto newparticle = world_.add_particle(product_species.id(), product_structure_id, new_com);
            if (rrec_) rrec_->StoreBindingReaction(time_, rrule.id(), reactant1.first, reactant2.first, newparticle.first);
            return std::vector<particle_id_pair>({ newparticle });
         }

         // Process (the lack of) change
         Log("GFRD").warn() << "pair reaction: placing product failed, moved reactants instead.";
         auto moved_particle1 = move_particle(reactant1, pos_pair1);
         auto moved_particle2 = move_particle(reactant2, pos_pair2);
         return std::vector<particle_id_pair>({ moved_particle1, moved_particle2 });
      }

      default:
         THROW_EXCEPTION(unsupported, "fire_pair_reaction: num products >= 2 not supported.");
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void process_single_event(Single& single, std::vector<DomainID>& ignore)
   {
      if (single.eventType() == Domain::EventType::INIT)             // SPECIAL CASE 1: If just need to make new domain.
      {
         make_new_domain(single);
      }
      else if (!std::isfinite(single.dt()))                          // SPECIAL CASE 2: In case nothing is scheduled to happen : do nothing and just reschedule
      {
         add_domain_event(single.id());            // this presumably only happens for non-defusing and non-decaying particles (D=0.0 and k=0.0)!
      }
      else                                                           // ### NORMAL: Process 'normal' event produced by the single
      {
         // 1. check that everything is ok
#if defined(_DEBUG)
         // The burst of the domain may be caused by an overlapping domain
         if (single.eventType() != Domain::EventType::BURST)
            ASSERT(check_domain(single));

         // check that the event time of the domain (last_time + dt) is equal to the simulator time
         ASSERT(std::abs(single.last_time() + single.dt() - time_) <= GfrdCfg.TIME_TOLERANCE * time_);
#endif

         // 2. get some necessary information
         auto pid_pair = single.pip();
         auto eventType = single.eventType();

         if (eventType == Domain::EventType::IV_EVENT)
            eventType = single.draw_iv_event_type(rng_);

         ASSERT(eventType == Domain::EventType::SINGLE_REACTION || eventType == Domain::EventType::SINGLE_ESCAPE || eventType == Domain::EventType::BURST);
         ReactionRule rrule = eventType == Domain::EventType::SINGLE_REACTION ? single.draw_reaction_rule(rng_) : ReactionRule(SpeciesTypeID(1), 0.0, std::vector<SpeciesTypeID>());      // ugly but only way to 

         // get the(new) position and structure on which the particle is located.
         position_structid_pair pos_struct = single.draw_new_position(rng_);
         pos_struct = world_.apply_boundary(pos_struct);

         // newpos now hold the new position of the particle (not yet committed to the world)
         // if we would here move the particles and make new shells, then it would be similar to a propagate

         remove_domain(single.id());      // single does not exist after this!!!

         std::vector<particle_id_pair> new_particles;
         switch (eventType) // possible actions for Single domains
         {
         case Domain::EventType::SINGLE_REACTION:
            new_particles = fire_single_reaction(pid_pair, rrule, pos_struct, ignore);
            break;

         case Domain::EventType::SINGLE_INTERACTION:
            //TODO particles, zero_singles_b, ignore = self.fire_interaction(single, newpos, struct_id, ignore)
            THROW_EXCEPTION(not_implemented, "SINGLE_INTERACTION not supported.")

         case Domain::EventType::SINGLE_ESCAPE:
         case Domain::EventType::BURST:
         {
            auto new_pid = fire_move(pid_pair, pos_struct);
            new_particles.emplace_back(new_pid);
         }  break;

         default:
            THROW_EXCEPTION(illegal_state, "unexpected eventType")
         }

         // 5. Make a new domain (or reuse the old one) for each particle(s)
         std::vector<Sphere> burst_spheres;
         for (auto& pid : new_particles)
         {
            auto did = create_single(pid);
            add_domain_event(did);
            // add to burst list
            burst_spheres.emplace_back(Sphere(pid.second.position(), pid.second.radius() * GfrdCfg.SINGLE_SHELL_FACTOR));              // Add the newly made init (zero-dt) singles to the list
            ignore.emplace_back(did);                       // Ignore these newly made singles(they should not be bursted)
         }

         // 6. Recursively burst around the newly made Init (zero-dt) NonInteractionSingles that surround the particles.
         // For each zero_single the burst radius equals the particle radius * SINGLE_SHELL_FACTOR
         // Add the resulting zero_singles to the already existing list.
         for (auto& s : burst_spheres)
         {
            burst_volume(s, ignore, true);
         }
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void process_pair_event(Pair& pair, std::vector<DomainID>& ignore)
   {
      // This method handles the things that need to be done when the current event was produced by a pair. The pair can be any type of pair(Simple or Mixed).
      // Note that this method is also called when a pair is bursted, in that case the event is a BURST event.
      // Note that ignore should already contain the id of 'pair'.

      // 1. check that everything is ok
#if defined(_DEBUG)
      // The burst of the domain may be caused by an overlapping domain
      if (pair.eventType() != Domain::EventType::BURST)
         ASSERT(check_domain(pair));

      // check that the event time of the domain (last_time + dt) is equal to the simulator time
      ASSERT(std::abs(pair.last_time() + pair.dt() - time_) <= GfrdCfg.TIME_TOLERANCE * time_);
#endif

      // 2. get some necessary information
      auto pid_pair1 = pair.pip1();
      auto pid_pair2 = pair.pip2();
      auto eventType = pair.eventType();

      // Draw actual pair event for iv at very last minute.
      if (eventType == Domain::EventType::IV_EVENT)
         eventType = pair.draw_iv_event_type(rng_);

      // 3.1 Get new position and current structures of particles
      position_structid_pair  pos_struct1, pos_struct2;
      if (pair.dt() > 0.0)
      {
         std::tie(pos_struct1, pos_struct2) = pair.draw_new_position(rng_, world_);
         pos_struct1 = world_.apply_boundary(pos_struct1);
         pos_struct2 = world_.apply_boundary(pos_struct2);

         //world_.check_overlap();
      }
      else
      {
         pos_struct1 = position_structid_pair(pid_pair1.second.position(), pid_pair1.second.structure_id());
         pos_struct2 = position_structid_pair(pid_pair2.second.position(), pid_pair2.second.structure_id());
      }

      // copy some more parameters from the pair before its destroyed.
      auto reactingsingle = pair.get_reacting_single();
      ReactionRule rrule = eventType == Domain::EventType::SINGLE_REACTION ? pair.draw_single_reaction_rule(rng_, reaction_rules_) :
         eventType == Domain::EventType::IV_REACTION ? pair.draw_reaction_rule(rng_) :
         ReactionRule(SpeciesTypeID(1), 0.0, std::vector<SpeciesTypeID>());      // ugly but only way to 
      Vector3 new_com = eventType == Domain::EventType::IV_REACTION ? world_.apply_boundary(pair.draw_new_com(rng_)) : Vector3();

      // newpos1 / 2 now hold the new positions of the particles(not yet committed to the world)
      // if we would here move the particles and make new shells, then it would be similar to a propagate

      remove_domain(pair.id());        // pair does not exist after this!!

      //3.2 If identity changing processes have taken place
      std::vector<particle_id_pair> new_particles;
      switch (eventType)
      {
      case Domain::EventType::SINGLE_REACTION:
         if (reactingsingle == pid_pair1.first)
         {
            auto ps = fire_move_pair(pid_pair2, pos_struct2, reactingsingle);
            new_particles = fire_single_reaction(pid_pair1, rrule, pos_struct1, ignore);
            new_particles.emplace_back(ps);
         }
         else
         {
            auto ps = fire_move_pair(pid_pair1, pos_struct1, reactingsingle);
            new_particles = fire_single_reaction(pid_pair2, rrule, pos_struct2, ignore);
            new_particles.emplace_back(ps);
         }
         break;

      case Domain::EventType::IV_REACTION:
         new_particles = fire_pair_reaction(pid_pair1, pos_struct1, pid_pair2, pos_struct2, rrule, new_com, ignore);
         break;

      case Domain::EventType::IV_ESCAPE:
      case Domain::EventType::COM_ESCAPE:
      case Domain::EventType::BURST:            // Just moving the particles
         new_particles.emplace_back(fire_move_pair(pid_pair1, pos_struct1, pid_pair2.first));
         new_particles.emplace_back(fire_move(pid_pair2, pos_struct2));
         break;

      default:
         THROW_EXCEPTION(illegal_state, "Unexpected event type");
      }

      // Put the zero singles resulting from putative bursting in fire_... in the final list.
      std::vector<Sphere> burst_spheres;
      for (auto& pid : new_particles)
      {
         auto did = create_single(pid);
         add_domain_event(did);
         // add to burst list
         burst_spheres.emplace_back(Sphere(pid.second.position(), pid.second.radius() * GfrdCfg.SINGLE_SHELL_FACTOR));              // Add the newly made init (zero-dt) singles to the list
         ignore.emplace_back(did);                       // Ignore these newly made singles(they should not be bursted)
      }

      // 6. Recursively burst around the newly made Init (zero-dt) NonInteractionSingles that surround the particles.
      // For each zero_single the burst radius equals the particle radius * SINGLE_SHELL_FACTOR
      // Add the resulting zero_singles to the already existing list.
      for (auto& s : burst_spheres)
      {
         burst_volume(s, ignore, true);
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void process_multi_event(Multi& multi, std::vector<DomainID>& ignore)
   {
      multi.step();
      switch (multi.eventType())
      {
      case Domain::EventType::MULTI_UNIMOLECULAR_REACTION:
      case Domain::EventType::MULTI_BIMOLECULAR_REACTION:
      case Domain::EventType::MULTI_ESCAPE:
         break_up_multi(multi, ignore);
         break;

      case Domain::EventType::MULTI_DIFFUSION:
         multi.set_last_time(time_);
         add_domain_event(multi.id());
         break;

      default: THROW_EXCEPTION(illegal_state, "Non multi EventType in multi event process");
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void break_up_multi(Multi& multi, std::vector<DomainID>& ignore)
   {
      // swap particles with multi
      std::vector<ParticleID> particles;
      particles.swap(const_cast<std::vector<ParticleID>&>(multi.get_particles()));

#if defined(_DEBUG)
      //Log("GFRD").warn(static_cast<std::string>(make_string() << "Multi existed from " << std::scientific << std::setprecision(9) << multi.start_time() << " to " << multi.last_time() << " approx. \t" << (uint)((multi.last_time() - multi.start_time())/ multi.dt()) << " steps, \t" << (uint)particles.size() << " particles.").c_str());
#endif

      // destroy domain
      remove_domain(multi.id());

      // create new init-single domains for each particle
      std::vector<Sphere> burst_spheres;
      for (auto pid : particles)
      {
         const auto& pid_pair = world_.get_particle(pid);
         auto did = create_single(pid_pair);
         add_domain_event(did);
         // add to burst list
         burst_spheres.emplace_back(Sphere(pid_pair.second.position(), pid_pair.second.radius() * GfrdCfg.SINGLE_SHELL_FACTOR));              // Add the newly made init (zero-dt) singles to the list
         ignore.emplace_back(did);                       // Ignore these newly made singles(they should not be bursted)
      }

      // 6. Recursively burst around the newly made Init (zero-dt) NonInteractionSingles that surround the particles.
      // For each zero_single the burst radius equals the particle radius * SINGLE_SHELL_FACTOR
      // Add the resulting zero_singles to the already existing list.
      for (auto& s : burst_spheres)
      {
         burst_volume(s, ignore, true);
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   // Recursively bursts the domains within 'radius' centered around 'pos'
   // Old code 'burst_non_multis' is replaced by this call, with skipMulties = true
   // Returns:
   // - the updated list domains that are already bursted -> ignore list
   // - zero_singles that were the result of the burst
   // NOTE: old used to multiply radius with SAFETY factor, do this now at the call-site
   void burst_volume(const Sphere& sphere, std::vector<DomainID>& ignore, bool skipMulties = false)
   {
      ShellCreateUtils::shell_burst_check<shell_matrix_type> sbc(sphere, ignore);
      CompileConfigSimulator::TBoundCondition::each_neighbor(shellmat_, sbc, sphere.position());
      auto to_burst = sbc.burst_range();

      for (auto did : to_burst)
      {
         if (std::find(ignore.cbegin(), ignore.cend(), did) != ignore.cend())          // if did to burst is on ignore list skip it (this could happen due to recursion!!) 
            continue;

         if (skipMulties)
         {
            auto& domain = domains_[did];
            if (domain->multiplicity() >= Domain::Multiplicity::MULTI) continue;    // skip multies
         }

         ignore.emplace_back(did);                                                     // put on ignore list

         // burst the domain, returns (one) init-single
         burst_domain(did, ignore);
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   // Reduces 'domain' (Single, Pair or Multi) to 'zero_singles', singles with the zero shell, and dt = 0 (INIT).
   // returns:
   // - list of zero_singles that was the result of the bursting
   // - updated ignore list
   void burst_domain(DomainID did, std::vector<DomainID>& ignore)
   {
      auto& domain = domains_[did];

      auto single = dynamic_cast<Single*>(domain.get());
      if (single != nullptr) return burst_single_domain(*single, ignore);

      auto pair = dynamic_cast<Pair*>(domain.get());
      if (pair != nullptr) return burst_pair_domain(*pair, ignore);

      auto multi = dynamic_cast<Multi*>(domain.get());
      if (multi != nullptr) return burst_multi_domain(*multi, ignore);

      THROW_EXCEPTION(not_implemented, "domain types not implemented.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void burst_single_domain(Single& single, std::vector<DomainID>& ignore)
   {
      // Check correct timeline 
      ASSERT(single.last_time() <= time_);
      ASSERT(time_ <= single.last_time() + single.dt());
      ASSERT(single.shell().code() == Shell::Code::NORMAL);       // cannot (or useless) burst an INIT shell

      // Override dt, the burst happens before the single's scheduled event.
      single.set_burst_time(time_);
      // with a burst there is always an associated event in the scheduler.
      // to simulate the natural occurence of the event we have to remove it from the scheduler
      scheduler_.remove(single.eventID());

      // calculate new_position, and make a INIT shell
      process_single_event(single, ignore);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void burst_pair_domain(Pair& pair, std::vector<DomainID>& ignore)
   {
      // Check correct timeline 
      ASSERT(pair.last_time() <= time_);
      ASSERT(time_ <= pair.last_time() + pair.dt());

      // Override dt, the burst happens before the single's scheduled event.
      pair.set_burst_time(time_);
      // with a burst there is always an associated event in the scheduler.
      // to simulate the natural occurence of the event we have to remove it from the scheduler
      scheduler_.remove(pair.eventID());

      // calculate new_position, and make a INIT shell
      process_pair_event(pair, ignore);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void burst_multi_domain(Multi& multi, std::vector<DomainID>& ignore)
   {
      // Check correct timeline 
      ASSERT(multi.last_time() <= time_);
      ASSERT(time_ <= multi.last_time() + multi.dt());

      // Override dt, the burst happens before the single's scheduled event.
      multi.set_burst_time(time_);
      // with a burst there is always an associated event in the scheduler.
      // to simulate the natural occurence of the event we have to remove it from the scheduler
      scheduler_.remove(multi.eventID());

      // calculate new_position, and make a INIT shells
      break_up_multi(multi, ignore);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   particle_id_pair fire_move(particle_id_pair pip, position_structid_pair pos_struct) const
   {
      // check if there is enough space
      auto ol = world_.test_particle_overlap(Sphere(pos_struct.first, pip.second.radius()), pip.first);     // check location of new position, ignore old location
      THROW_UNLESS_MSG(no_space, !ol, "No space to move particle:" << pip.first);
      return move_particle(pip, pos_struct);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   particle_id_pair fire_move_pair(particle_id_pair pip, position_structid_pair pos_struct, const ParticleID partner_id) const
   {
      // check if there is enough space
      auto ol = world_.test_particle_overlap(Sphere(pos_struct.first, pip.second.radius()), std::vector<ParticleID>{ pip.first, partner_id });     // check location of new position, ignore old location and partner particle
      THROW_UNLESS_MSG(no_space, !ol, "No space to move particle:" << pip.first << " partner:" << partner_id);
      return move_particle(pip, pos_struct);
   };

   // --------------------------------------------------------------------------------------------------------------------------------

   particle_id_pair move_particle(particle_id_pair pip, position_structid_pair pos_struct) const
   {
      // Moves a particle in World and performs the required change of structure_id based on an existing particle.
      auto& p = pip.second;
      particle_id_pair new_pid_particle_pair = std::make_pair(pip.first, Particle(p.sid(), Sphere(pos_struct.first, p.radius()), p.structure_id(), p.D(), p.v()));
      world_.update_particle(new_pid_particle_pair);
      return new_pid_particle_pair;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void remove_domain(DomainID did)
   {
      // Removes all the ties to a domain(single, pair, multi) from the system.
      // Note that the particles that it represented still exist in 'world'
      auto& domain = domains_[did];
      if (domain->num_shells() == 1)
      {
         auto shell_id_ref = domain->get_shell();
         shellmat_.erase(shell_id_ref.first);
      }
      else
      {
         for (auto shell_id_ref : domain->get_shell_list())
            shellmat_.erase(shell_id_ref.first);
      }

      domains_.erase(did);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void make_new_domain(Single& single)
   {
      THROW_UNLESS(illegal_state, single.eventType() == Domain::EventType::INIT);

      auto single_pos = single.particle().position();

      // collect all distances and radii of surrounding particles with init-shells

      bool succes;
      ShellCreateUtils::shell_interaction_check<shell_matrix_type> sic(single.shell_id(), single.particle());
      CompileConfigSimulator::TBoundCondition::each_neighbor(shellmat_, sic, single_pos);
      if (sic.multiple())
      {
         // make a MULTI domain
         succes = form_multi(single.id(), sic.multi_range());
      }
      else if (sic.did())
      {
         // Make a PAIR domain, with this single and sic.did()
         succes = try_pair(single.id(), sic.did());

         // When pair failed, try create a Single
         if (!succes)
         {
            succes = update_single(single);
         }
      }
      else
      {
         // Just scale and update the Single domain
         succes = update_single(single);
      }

      if (!succes)
      {
         succes = form_multi(single.id(), sic.multi_range());
         THROW_UNLESS_MSG(illegal_state, succes, "Failed to make new domain!");
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool try_pair(DomainID did1, DomainID did2)
   {
      // Try to make a pair domain out of the two singles.A pair domain can be a :
      // -SimplePair (both particles live on the same structure)
      // -MixedPair (the particles live on different structures -> MixedPair2D3D)
      // Note that single1 is always the single that initiated the creation of the pair
      //  and is therefore no longer in the scheduler

      Single* single1 = dynamic_cast<Single*>(domains_[did1].get());
      Single* single2 = dynamic_cast<Single*>(domains_[did2].get());
      particle_id_pair pip1 = single1->pip();
      particle_id_pair pip2 = single2->pip();

      DomainID did = didgen_();
      ShellID sid = sidgen_();         // TODO, maybe move id-gen to after we know there is space!

      const auto& rules = reaction_rules_.query_reaction_rules(pip1.second.sid(), pip2.second.sid());
      auto structure1 = world_.get_structure(pip1.second.structure_id());
      auto structure2 = world_.get_structure(pip2.second.structure_id());
      auto pair = create_default_pair(did, sid, pip1, structure1, pip2, structure2, rules);
      if (pair != nullptr)
      {
         bool succes = pair->create_updated_shell(shellmat_, world_, single1->shell_id(), single2->shell_id());
         if (!succes) return false;

         pair->set_k_totals(single1->particle().sid(), single1->k_total(), single2->k_total());
         // # 3. update the shell containers with the new shell
         shellmat_.update(pair->shell_pair());

         // determine next action for this event
         pair->determine_next_event(rng_);
         pair->set_last_time(time_);

         // clean old single2 event
         scheduler_.remove(single2->eventID());

         // remove singles and add new pair domain
         remove_domain(did1);
         remove_domain(did2);
         domains_[did] = std::unique_ptr<Domain>(std::move(pair));

         // update scheduler
         add_domain_event(did);
      }
      else
      {
         Log("GFRD").info() << "Make Pair from " << did1 << " and " << did2 << " failed (combination of surfaces not implemented)!";
         return false;
      }
      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool form_multi(DomainID did, const gi::iteratorRange<std::vector<DomainID>>& neighbors)
   {
      Multi* multi = nullptr;
      std::vector<DomainID> neighbor_list;

      // Copy known neigbors to our neigbor list
      std::copy(neighbors.begin(), neighbors.end(), std::back_inserter(neighbor_list));

      // loop all neighbor domains, add any init/multi domains whos particles lays within the multi_horizon distance
      {
         std::vector<DomainID> check_neighbors;
         std::copy(neighbor_list.begin(), neighbor_list.end(), std::back_inserter(check_neighbors));
         while (!check_neighbors.empty())
         {
            std::vector<DomainID> more_neighbors;
            for (auto ndid : check_neighbors)
            {
               const auto& domain = domains_[ndid];

               // is new_neighbor a single(init)?
               auto single = dynamic_cast<Single*>(domain.get());
               if (single != nullptr)
               {
                  // make room for MSF shell ?        TODO not sure if this is needed, and its also expensive since we also do sic here.
                  // ignore.emplace_back(ndid);
                  // burst_volume(Sphere(single->particle().position(), single->particle().radius() * GfrdCfg.MULTI_SHELL_FACTOR * GfrdCfg.SAFETY), ignore, true);

                  // add other (init)singles and multies within multi_horizon, that are not yet in the list
                  ShellCreateUtils::shell_interaction_check<shell_matrix_type> sic(single->shell_id(), single->particle());
                  CompileConfigSimulator::TBoundCondition::each_neighbor(shellmat_, sic, single->particle().position());
                  for (auto edid : sic.multi_range())
                  {
                     // already in list?
                     if (did != edid &&
                        std::end(neighbor_list) == std::find(neighbor_list.begin(), neighbor_list.end(), edid) &&
                        std::end(check_neighbors) == std::find(check_neighbors.begin(), check_neighbors.end(), edid) &&
                        std::end(more_neighbors) == std::find(more_neighbors.begin(), more_neighbors.end(), edid))
                        more_neighbors.emplace_back(edid);
                  }
                  continue;
               }

               // is new_neighbor a multi?
               auto multi1 = dynamic_cast<Multi*>(domain.get());
               if (multi1 != nullptr)
               {
                  if (multi == nullptr || multi1 == multi)      // first multi in list (reuse it), or re-hit the same (remove from list)
                  {
                     multi = multi1;         // use exiting multi, remove from neighbor list
                     neighbor_list.erase(std::find(neighbor_list.begin(), neighbor_list.end(), ndid));
                  }
               }
            }

            // add new particles to neighbor list
            std::copy(more_neighbors.begin(), more_neighbors.end(), std::back_inserter(neighbor_list));
            // next loop over new particles
            check_neighbors.clear();
            std::copy(more_neighbors.begin(), more_neighbors.end(), std::back_inserter(check_neighbors));
         }
      }

      // create a new multi if needed
      bool new_multi = false;
      if (multi == nullptr)
      {
         new_multi = true;
         DomainID newdid = didgen_();
         auto domain = std::make_unique<Multi>(newdid, world_, shellmat_, reaction_rules_, rng_, rrec_, time_);
         multi = domain.get();
         domains_[newdid] = std::unique_ptr<Domain>(std::move(domain));
      }

      // put original Single into Multi
      auto single = dynamic_cast<Single*>(domains_[did].get());
      add_to_multi(multi, single->pip());
      remove_domain(did);

      // Add our neighbor particles, and merge other multi's.
      for (auto ndid : neighbor_list)
      {
         auto multi1 = dynamic_cast<Multi*>(domains_[ndid].get());
         if (multi1 != nullptr)
         {
            THROW_UNLESS_MSG(illegal_state, multi1 != multi, "Cannot add multi ourself to ourself!");
            // copy particles and shells from other multi to this multi (don't move the shells, it may cause overlap)
            add_multi_to_multi(multi, multi1);

            scheduler_.remove(multi1->eventID());
            remove_domain(ndid);
            continue;
         }

         auto single1 = dynamic_cast<Single*>(domains_[ndid].get());
         if (single1 != nullptr)
         {
            add_to_multi(multi, single1->pip());
            // clean old single event and remove domain
            scheduler_.remove(single1->eventID());
            remove_domain(single1->id());
         }
      }

      multi->initialize(time());
      if (!new_multi)
         scheduler_.update(multi->eventID(), Event(time_ + multi->dt(), multi->id()));         // reschedule
      else
         add_domain_event(multi->id());                                                        // add to scheduler

      // TODO check why so many single particle multies are created
      //size_t count = multi->get_particles().size();
      //if (count < 2) dump(make_string() << "d:\\small_multi.log", true);

      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void add_to_multi(Multi* multi, particle_id_pair pid_pair)
   {
      multi->queue_particle(pid_pair);

      const auto sid_pair = std::make_pair<ShellID, Shell>(sidgen_(), Shell(multi->id(), Sphere(pid_pair.second.position(), pid_pair.second.radius() * GfrdCfg.MULTI_SHELL_FACTOR), Shell::Code::MULTI));
      shellmat_.update(sid_pair);
      multi->add_shell(sid_pair, pid_pair.first);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void add_multi_to_multi(Multi* multi, Multi* source)
   {
      for (auto pid : source->get_particles())
      {
         const auto& pip = world_.get_particle(pid);
         multi->queue_particle(pip);

         const auto& sid_old = source->get_shell_for_particle(pip.first);     // new shell but in old location
         const auto sid_pair = std::make_pair<ShellID, Shell>(sidgen_(), Shell(multi->id(), sid_old.get_sphere(), Shell::Code::MULTI));
         shellmat_.update(sid_pair);
         multi->add_shell(sid_pair, pip.first);
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool update_single(Single& single)
   {
      if (!single.create_updated_shell(shellmat_, world_)) return false;

      shellmat_.update(single.shell_pair());
      single.determine_next_event(rng_);
      single.set_last_time(time_);

      add_domain_event(single.id());
      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void check_intern() const
   {
      world_.check();

      THROW_UNLESS_MSG(illegal_state, scheduler_.check(), "Schedular check failed.");
      THROW_UNLESS_MSG(illegal_state, time_ >= 0.0, "simulation time negative (t=" << time_ << ").");
      THROW_UNLESS_MSG(illegal_state, dt_ >= 0.0, "simulation timestep negative (dt=" << dt_ << ").");

      size_t shellcount = 0;
      for (const auto& domain : domains_)
      {
         THROW_UNLESS_MSG(illegal_state, domain.second, "DomainID:" << domain.first << " ptr is null.");
         shellcount += check_domain(*domain.second.get());

         auto eventID = domain.second.get()->eventID();
         THROW_UNLESS_MSG(illegal_state, scheduler_.has(eventID), "Domain Event is not in the schedular.");
         auto ev = scheduler_.get(eventID);      // domain has no event
         THROW_UNLESS_MSG(illegal_state, ev.time() >= time_, "Domain Event time should lay in the future.");
      }

      // check particle positions and overlap
      check_particles();

      // some checks are valid only after first step
      if (num_steps_ == 0) return;

      // number of shells in domains should equal total shells in shell_matrix
      THROW_UNLESS_MSG(illegal_state, shellmat_.size() == shellcount, "Shells in matrix " << shellmat_.size() << " does not match shell count from domains " << shellcount << ".");

      // count shells (pair=2) and compare to number of particles
      size_t pcount = std::accumulate(domains_.begin(), domains_.end(), size_t(0), [](size_t sum, const domain_map::value_type& d) { return sum + (d.second->multiplicity() == Domain::Multiplicity::PAIR ? 2 : d.second->num_shells()); });
      THROW_UNLESS_MSG(illegal_state, pcount == world_.num_particles(), "Counted shells " << pcount << " does not match number of particles " << world_.num_particles() << ".");


      // protective domains (shells) cannot overlap, except multies
      check_shells();

      // TODO
      /*

      self.check_shell_matrix()
      self.check_domains()
      self.check_event_stoichiometry()
      self.check_particle_consistency()

      */

   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void check_shells() const
   {
      for (const auto& i : shellmat_)
      {
         const auto& shell = i.second;
         DomainID did = shell.did();
         auto& domain = domains_.at(did);
         std::vector<DomainID> ovl;

         if (shell.shape() == Shell::Shape::SPHERE)
         {
            auto sphere = shell.get_sphere();
            ShellCreateUtils::shell_overlap_check<shell_matrix_type> soc(sphere);
            CompileConfigSimulator::TBoundCondition::each_neighbor(shellmat_, soc, sphere.position());
            ovl = std::move(soc.overlap());
         }
         else
         {
            //auto cylinder = shell.get_cylinder();
            //ShellCreateUtils::shell_overlap_check<shell_matrix_type> soc(cylinder);       // TODO check cylinder overlaps
            //CompileConfigSimulator::TBoundCondition::each_neighbor(shellmat_, soc, cylinder.position());
            //ovl = std::move(soc.overlap());
         }

         auto multi = dynamic_cast<Multi*>(domain.get());
         if (multi != nullptr)
         {
            THROW_UNLESS_MSG(illegal_state, ovl.size() >= 1 && ovl.size() <= multi->num_shells(), "Shells of multi domain:" << did << " overlaps with more shells (" << ovl.size() << ") than it holds (" << multi->num_shells() << ").");
            for (const auto& j : ovl)
               THROW_UNLESS_MSG(illegal_state, j == did, "Shell of multi domain:" << did << " overlaps with shell from domain:" << j);
         }
         else
         {
            THROW_UNLESS_MSG(illegal_state, ovl.size() == 1, "Shell of domain:" << did << " overlaps with other shells.");
            THROW_UNLESS_MSG(illegal_state, *ovl.begin() == did, "Shell of domain:" << did << " overlaps with other shell (Domain:" << *ovl.begin() << " ).");
         }
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void check_particles() const
   {
      // check overlaping (or touching) particles at start (user input error)
      for (const auto& p : world_.get_particles())
      {
         Vector3 ws = world_.world_size();
         Vector3 pos = p.second.position();
         THROW_UNLESS_MSG(illegal_state, pos.X() >= 0 && pos.X() < ws.X(), "Particle:" << p.first << " is out of this world! X:" << pos.X() << " WSx: 0 - " << ws.X() << ".");
         THROW_UNLESS_MSG(illegal_state, pos.X() >= 0 && pos.X() < ws.X(), "Particle:" << p.first << " is out of this world! Y:" << pos.Y() << " WSy: 0 - " << ws.Y() << ".");
         THROW_UNLESS_MSG(illegal_state, pos.Z() >= 0 && pos.Z() < ws.Z(), "Particle:" << p.first << " is out of this world! Z:" << pos.Z() << " WSz: 0 - " << ws.Z() << ".");
         THROW_UNLESS_MSG(illegal_state, !world_.test_particle_overlap(p.second.shape(), p.first), "Particle:" << p.first << " overlaps with an other particle.");
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   size_t check_domain(const Domain& domain) const         // may need to split into Single/Pair/Multi checks
   {
      auto single = dynamic_cast<const Single*>(&domain);
      if (single != nullptr)
      {
         auto shell = single->shell();
         check_shell(shell);
         return 1;
      }

      auto pair = dynamic_cast<const Pair*>(&domain);
      if (pair != nullptr)
      {
         auto shell = pair->shell();
         check_shell(shell);
         return 1;
      }

      auto multi = dynamic_cast<const Multi*>(&domain);
      if (multi != nullptr)
      {
         for (auto& sid_pair : multi->get_shell_list())
         {
            auto& shell = sid_pair.second.get();
            THROW_UNLESS_MSG(illegal_state, shell.shape() == Shell::Shape::SPHERE, "Domain:" << domain.id() << " is a multi but not Spherical.");        // only spherical shells in a multi
            check_shell(shell);
         }
         return multi->get_shell_list().size();
      }

      // TODO egfrd.py:3093
      return 0;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void check_shell(const Shell& shell) const
   {
      World::particle_id_pair_and_distance_list ovl;

      if (shell.shape() == Shell::Shape::SPHERE)
      {
         auto s = shell.get_sphere();
         THROW_UNLESS_MSG(illegal_state, s.radius() <= shellmat_.cell_size() / 2.0, "Shell (Sphere) radius may not exceed half the matrix cell size.");

         // find particles within this protective sphere
         ovl = std::move(world_.check_particle_overlap(s));
      }
      else
      {
         auto c = shell.get_cylinder();
         THROW_UNLESS_MSG(illegal_state, Vector2(c.half_length(), c.radius()).length() < shellmat_.cell_size() / 2.0, "Shell (Cylinder) any angle (corner distance) may not exceed half the matrix cell size.");

         // find particles within this protective cylinder
         ovl = std::move(world_.check_particle_overlap(c));
      }

      DomainID did = shell.did();
      auto& domain = domains_.at(did);

      // SINGLE
      auto single = dynamic_cast<Single*>(domain.get());
      if (single != nullptr)
      {
         THROW_UNLESS_MSG(illegal_state, ovl.size() == 1, "Single Domain:" << did << " shell may enclose only one particle.");
         ParticleID p = (*ovl.begin()).first.first;
         THROW_UNLESS_MSG(illegal_state, ovl[0].first.first == single->particle_id(), "Single Domain:" << did << " shell should enclose Particle:" << single->particle_id() << " but it encloses Particle:" << p << ".");
         if (shell.code() == Shell::Code::INIT)
         {
            THROW_UNLESS_MSG(illegal_state, -ovl[0].second == 2 * single->particle().radius(), "Single Init Domain:" << did << " shell size does not match Particle:" << single->particle_id() << " size.");
         }
         else
         {
            THROW_UNLESS_MSG(illegal_state, -ovl[0].second > 2 * single->particle().radius(), "Single Domain:" << did << " has Particle:" << single->particle_id() << " sticking out of its shell.");
         }
      }

      // PAIR
      auto pair = dynamic_cast<Pair*>(domain.get());
      if (pair != nullptr)
      {
         THROW_UNLESS_MSG(illegal_state, ovl.size() == 2, "Pair Domain:" << did << " shell may enclose only two particles.");
         ParticleID p1 = ovl[0].first.first;
         ParticleID p2 = ovl[1].first.first;
         THROW_UNLESS_MSG(illegal_state, (p1 == pair->particle1_id() && p2 == pair->particle2_id()) || (p2 == pair->particle1_id() && p1 == pair->particle2_id()),
            "Pair Domain:" << did << " shell should enclose Particles:" << pair->particle1_id() << " & " << pair->particle2_id() << " but it encloses Particles:" << p1 << " & " << p2 << ".");
         if (p1 == pair->particle1_id())
         {
            THROW_UNLESS_MSG(illegal_state, -ovl[0].second > 2 * pair->particle1().radius(), "Pair Domain:" << did << " has Particle:" << pair->particle1_id() << " sticking out of its shell.");
            THROW_UNLESS_MSG(illegal_state, -ovl[1].second > 2 * pair->particle2().radius(), "Pair Domain:" << did << " has Particle:" << pair->particle2_id() << " sticking out of its shell.");
         }
         else
         {
            THROW_UNLESS_MSG(illegal_state, -ovl[0].second > 2 * pair->particle2().radius(), "Pair Domain:" << did << " has Particle:" << pair->particle2_id() << " sticking out of its shell.");
            THROW_UNLESS_MSG(illegal_state, -ovl[1].second > 2 * pair->particle1().radius(), "Pair Domain:" << did << " has Particle:" << pair->particle1_id() << " sticking out of its shell.");
         }
      }

      // MULTI
      auto multi = dynamic_cast<Multi*>(domain.get());
      if (multi != nullptr)
      {
         const auto& particles = multi->get_particles();
         THROW_UNLESS_MSG(illegal_state, ovl.size() > 0 && ovl.size() <= multi->num_shells(), "Multi Domain:" << did << " shell may enclose between 1 and " << multi->num_shells() << " particles.");
         int shellcount = 0;
         for (auto& i : ovl)
         {
            ParticleID p = i.first.first;
            const auto& fp = std::find(particles.begin(), particles.end(), p);
            THROW_UNLESS_MSG(illegal_state, fp != particles.end(), "Multi Domain:" << did << " shell encloses Particle:" << p << " that does not belong to the multi.");
            if (multi->get_shell_for_particle(p) == shell)
            {
               shellcount++;
               const auto& pp = world_.get_particle(p);
               THROW_UNLESS_MSG(illegal_state, -ovl[0].second > 2 * pp.second.radius(), "Multi Domain:" << did << " has Particle:" << p << " sticking out of its shell.");
            }
         }
         THROW_UNLESS_MSG(illegal_state, shellcount == 1, "Multi Domain:" << did << " shell links to multiple particles.");      // link this shell to only one particle in multi
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   EGFRDSimulator() = delete;
   EGFRDSimulator(const EGFRDSimulator&) = delete;                                     // don't allow value construct of the simulator
   EGFRDSimulator& operator=(const EGFRDSimulator& w) = delete;                        // don't allow value assignment of the simulator
   EGFRDSimulator(EGFRDSimulator&&) = default;                                         // allow move construction of a simulator (note that it implementation may be a copy!)
   EGFRDSimulator& operator=(EGFRDSimulator&& w) = default;                            // allow move assignment of a simulator

   // --------------------------------------------------------------------------------------------------------------------------------

protected:
   friend class Persistence;

   SerialIDGenerator<ShellID> sidgen_;
   SerialIDGenerator<DomainID> didgen_;
   EventScheduler scheduler_;

   domain_map domains_;
   shell_matrix_type shellmat_;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* EGFRD_SIMULATOR_HPP */