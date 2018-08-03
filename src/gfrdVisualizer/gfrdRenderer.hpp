#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <World.hpp>
#include <EGFRDSimulator.hpp>
#include "CameraController.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class ExtSim;

// --------------------------------------------------------------------------------------------------------------------------------

inline void glVertex3(Vector3 v) { glVertex3d(v.X(), v.Y(), v.Z()); }

void RenderGFRD(const EGFRDSimulator& s, const CameraController& cam, int showpid, bool drawShells, StructureTypeID selStructure, ParticleID selParticle, DomainID selDomain, std::array<uint, 3>& domain_type_count);
void RenderGFRDExtern(const ExtSim& w, const CameraController& cam, int showpid, bool drawShells, uint selStructure, uint selParticle, std::array<uint, 3>& domain_type_count);

void drawAllStructures(StructureContainer::structures_range, const CameraController& cam, StructureTypeID selStructure);
void drawAllShells(const std::vector<std::pair<std::reference_wrapper<const Shell>, uint>>& shells, const Vector3& u);
void drawMatrix(double cs, std::array<uint, 3> ms);
void drawParticle(const Particle::particle_id_pair& p, ParticleID selParticle, int showID, const Vector3& v, const World* world);

void drawBox(const Box& box);
void drawPlane(const Plane& plane, const StructureTypeID& id, bool select);
void drawDisk(const Disk& disk, const StructureTypeID& id, bool select);
void drawSphere(const Sphere& shpere, const StructureTypeID& sid, bool select);
void drawCylinder(const Cylinder& cylinder, const StructureTypeID& sid, bool select);

// --------------------------------------------------------------------------------------------------------------------------------

void glxSolidCylinder(double radius, double length, int sides);
void glxWireCylinder(double radius, double length, int sides);

// --------------------------------------------------------------------------------------------------------------------------------

