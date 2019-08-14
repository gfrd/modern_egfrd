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
    }

    ~StructureSection() = default;

    // --------------------------------------------------------------------------------------------------------------------------------

    static std::string section_name() { return "Structure"; }

    const std::string key_name = "Name";
    const std::string &name() const { return name_; }

    const std::string key_type = "Type";
    const std::string& structureType() const { return structure_type_; }

    const std::string key_pos = "Pos";
    Vector3 pos() const { return pos_ }


    // --------------------------------------------------------------------------------------------------------------------------------

    bool set_keypair(const std::string &key, const std::string &value) override
    {
        if (SectionBase::set_keypair(key, value)) return true;
        if (key == key_name)
        {
            THROW_UNLESS_MSG(illegal_section_value, name_.empty(),
                             "Name already set! Use a new Structure section to define a new structure.");
            name_ = format_check(value);
            return true;
        }
        if (key == key_pos)
        {
            bool found;
            pos_ = get_vec3(value, "Structure position", (bool &) &found);
            return true;
        }
        THROW_EXCEPTION(illegal_section_key, "Key '" << key << "' not recognized.");
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    void create_structure(Model &model, const VariablesSection &vars) const
    {
        auto sid = model.get_def_structure_type_id();
        for (auto name : names_)
        {
//            model.add_species_type(SpeciesType(name, sid, D(), r(), v()));
        }
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

    std::string format_check(std::string input)
    {
        auto regex = std::regex("[^\\s,]+");
        auto begin = std::sregex_iterator(input.begin(), input.end(), regex);
        auto end = std::sregex_iterator();

        auto size = std::distance(begin, end);
        THROW_UNLESS_MSG(illegal_section_value, size > 0, "SpeciesTypeName not valid");
        names_.reserve(size);

        std::stringstream ss;
        for (auto i = begin; i != end; ++i)
        {
            std::string name = i->str();
            if (!is_valid_speciestype_name(name)) THROW_EXCEPTION(illegal_section_value,
                                                                  "SpeciesTypeName '" << name << "'not valid");
            names_.emplace_back(name);
            if (i != begin) ss << ", ";
            ss << name;
        }

        return ss.str();
    }

    // --------------------------------------------------------------------------------------------------------------------------------

    std::string name_, structure_type_;
    std::vector<std::string> names_;
    Vector3 pos_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream &operator<<(std::ostream &stream, const StructureSection &sts)
{
    stream << "[" << sts.section_name() << "]" << std::endl;
    stream << sts.key_name << " = " << sts.name() << std::endl;
    stream << sts.key_radius << " = " << sts.r() << std::endl;
    stream << sts.key_diffusion << " = " << sts.D() << std::endl;
    if (sts.v() != 0) stream << sts.key_drift_velocity << " = " << sts.v() << std::endl;
    stream << sts.key_structure_type << " = " << sts.structureType() << std::endl;
    stream << std::endl;
    return stream;
}
