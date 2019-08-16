#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "SpeciesType.hpp"
#include "Model.hpp"
#include "SectionBase.hpp"
#include <iomanip>
#include <Vector3.hpp>

// --------------------------------------------------------------------------------------------------------------------------------

struct ME_EXPORT StructureSection final : SectionBase
{
    explicit StructureSection() : SectionBase()
    {
        init_auto_vars( { { key_lx, 0.0}, { key_ly, 0.0} } );
    }

    ~StructureSection() = default;

    // --------------------------------------------------------------------------------------------------------------------------------

    static std::string section_name() { return "Structure"; }

    const std::string key_name = "Name";
    const std::string &name() const { return name_; }

    const std::string key_type = "Type";
    const std::string& structureType() const { return structure_type_; }

    const std::string key_pos = "Position";
    Vector3 pos() const { return pos_; }

    const std::string key_vx = "Vx";
    Vector3 vx() const { return vx_; }

    const std::string key_vy = "Vy";
    Vector3 vy() const { return vy_; }

    const std::string key_lx = "Lx";
    double lx() const { return auto_var_value(key_lx); }

    const std::string key_ly = "Ly";
    double ly() const { return auto_var_value(key_ly); }


    // --------------------------------------------------------------------------------------------------------------------------------

    bool set_keypair(const std::string &key, const std::string &value) override
    {
        if (SectionBase::set_keypair(key, value)) return true;
        if (key == key_name)
        {
            THROW_UNLESS_MSG(illegal_section_value, name_.empty(),
                             "Name already set! Use a new Structure section to define a new structure.");
            name_ = value;
            return true;
        }
        if (key == key_type)
        {
            THROW_UNLESS_MSG(illegal_section_value, structure_type_.empty(),
                             "Type already set! Use a new Structure section to define a new structure.");
            structure_type_ = value;
            return true;
        }
        if (key == key_pos)
        {
            bool found;
            pos_ = get_vec3(value, "Structure position", found);

            THROW_UNLESS_MSG(illegal_section_value, found,
                             "Structure position was not recognised. Ensure the format is (x, y, z), where each coordinate can be a value or an expression.");
            return true;
        }
        if (key == key_vx)
        {
            bool found;
            vx_ = get_vec3(value, "Structure Vx", found);

            THROW_UNLESS_MSG(illegal_section_value, found,
                             "Structure Vx was not recognised. Ensure the format is (x, y, z), where each coordinate can be a value or an expression.");
            return true;
        }
        if (key == key_vy)
        {
            bool found;
            vy_ = get_vec3(value, "Structure Vy", found);

            THROW_UNLESS_MSG(illegal_section_value, found,
                             "Structure Vy was not recognised. Ensure the format is (x, y, z), where each coordinate can be a value or an expression.");
            return true;
        }
        THROW_EXCEPTION(illegal_section_key, "Key '" << key << "' not recognized.");
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    void create_structure_type(Model &model) const
    {
        if(model.get_structure_type_id_by_name(name_) == StructureTypeID(0))
        {
            // Only add structure type if it doesn't exist yet
            model.add_structure_type(StructureType(name_));
        }
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    void create_structure(Model &model, World &world) const
    {
        auto wsid = world.get_def_structure_id();
        auto sid = model.get_structure_type_id_by_name(name_);

        THROW_UNLESS_MSG(illegal_section_value, sid != StructureTypeID(0), "Structure type '" << structure_type_ << "' not recognized.")

        if (structure_type_ == "PlanarSurface") {
            auto plane = Plane(pos_, vx_, vy_, lx(), ly(), false);
            auto structure = std::make_shared<PlanarSurface>(PlanarSurface(name_, sid, wsid, plane));
            world.add_structure(structure);
            return;
        }

        THROW_EXCEPTION(illegal_section_value, "Structure type '" << structure_type_ << "' not recognized.");
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    Vector3 get_vec3(const std::string& value, const std::string name, bool& found)
    {
        std::smatch match;
        auto regex = std::regex(R"(\(([^,]+),([^,]+),([^,]+)\))");
        found = std::regex_search(value, match, regex);
        if (!found || match.size() != 4) {
            found = false;
            return Vector3();
        }

        auto x = vars_->evaluate_value_expression(match[1].str(), name);
        auto y = vars_->evaluate_value_expression(match[2].str(), name);
        auto z = vars_->evaluate_value_expression(match[3].str(), name);

        return Vector3(x, y, z);
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    void PrintSettings() const override
    {
//        std::cout << std::setw(14) << "species = " << "'" << name_ << "'" << ", D = " << D() << " [m^2*s^-1], r = "
//                  << r() << " [m]";
//        if (v() != 0) std::cout << ", v = " << v();
//        std::cout << "\n";
    }

    // --------------------------------------------------------------------------------------------------------------------------------

private:

    // --------------------------------------------------------------------------------------------------------------------------------

    std::string name_, structure_type_;
    Vector3 pos_, vx_, vy_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream &operator<<(std::ostream &stream, const StructureSection &sts)
{
    stream << "[" << sts.section_name() << "]" << std::endl;
    stream << sts.key_name << " = " << sts.name() << std::endl;
//    stream << sts.key_radius << " = " << sts.r() << std::endl;
//    stream << sts.key_diffusion << " = " << sts.D() << std::endl;
//    if (sts.v() != 0) stream << sts.key_drift_velocity << " = " << sts.v() << std::endl;
//    stream << sts.key_structure_type << " = " << sts.structureType() << std::endl;
    stream << std::endl;
    return stream;
}
