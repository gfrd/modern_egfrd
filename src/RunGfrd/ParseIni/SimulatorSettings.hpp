#ifndef SIMULATORSETTINGS_HPP
#define SIMULATORSETTINGS_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <vector>
#include "SimulatorSection.hpp"
#include "WorldSection.hpp"
#include "SpeciesTypeSection.hpp"
#include "ReactionRuleSection.hpp"
#include "ParticlesSection.hpp"
#include "ParserExceptions.hpp"
#include "CopyNumbersSection.hpp"
#include "ParticlePositionsSection.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class SimulatorSettings
{
public:
   explicit SimulatorSettings() : current_(nullptr) { }
   ~SimulatorSettings() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   void create_section_heading(const std::string line)
   {
      std::smatch match;
      bool is_match = std::regex_search(line, match, regex_pattern_);
      THROW_UNLESS_MSG(illegal_section, is_match, line);

      auto heading = match[1].str();
      create_section(heading);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool is_section_heading(const std::string& line) const
   {
      return std::regex_match(line, regex_pattern_);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void add_keypair(const std::string& line) const
   {
      keypair_from_line(line);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:

   void keypair_from_line(const std::string& line) const
   {
      auto key = std::string();
      auto value = std::string();
      bool is_match = get_name_value(line, key, value);

      if (is_match && current_ != nullptr) current_->set_keypair(key, check_parameter(value));
   }

   std::string check_parameter(std::string value) const
   {
      std::smatch match;
      auto regex = std::regex("<param(\\d)>");
      bool is_match = std::regex_search(value, match, regex);
      if (!is_match) return value;
      auto index = match[1].str();
      int p = std::stoi(index);
      return parameters_[p];
   }

   bool get_name_value(const std::string& line, std::string& name, std::string& value) const
   {
      std::smatch match;
      auto regex = std::regex("^\\s*(\\w+)\\s*=\\s*([\\w<>]+[^;]*).*");
      bool is_match = std::regex_search(line, match, regex);
      name = match[1].str();
      value = match[2].str();
      return is_match;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void create_section(const std::string& heading)
   {
      if (heading == SimulatorSection::section_name()) { current_ = &simulatorSection_; return; }
      if (heading == WorldSection::section_name()) { current_ = &worldSection_; return; }
      if (heading == ParticlesSection::section_name()) { current_ = &throwInParticlesSection_; return; }
      if (heading == CopyNumbersSection::section_name()) { current_ = &copyNumbersSection_; return; }
      if (heading == ParticlePositionSection::section_name()) { current_ = &particlePositionsSection_; return; }
      if (heading == ReactionRuleSection::section_name()) { reactionRuleSections_.emplace_back(); current_ = &reactionRuleSections_[reactionRuleSections_.size() - 1]; return; }
      if (heading == SpeciesTypeSection::section_name()) { speciesTypeSections_.emplace_back(); current_ = &speciesTypeSections_[speciesTypeSections_.size() - 1]; return; }
      THROW_EXCEPTION(illegal_section, "Section '" << heading << "' not recognized!");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

public:

   const SimulatorSection& getSimulatorSection() const { return simulatorSection_; }
   const WorldSection& getWorldSection() const { return worldSection_; }
   const std::vector<SpeciesTypeSection>& getSpeciesTypeSections() const { return speciesTypeSections_; }
   const std::vector<ReactionRuleSection>& getReactionRuleSections() const { return reactionRuleSections_; }
   const ParticlesSection& getParticlesSection() const { return throwInParticlesSection_; }
   const CopyNumbersSection& getCopyNumbersSection() const { return copyNumbersSection_; }
   const ParticlePositionSection& getParticlePositionsSection() const { return particlePositionsSection_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   size_t parameter_size() const { return parameters_.size(); }
   void set_parameter(size_t i, std::string value) { parameters_[i] = value; }
   const std::string& get_parameter(size_t i) { return parameters_[i]; }

   // --------------------------------------------------------------------------------------------------------------------------------

private:

   SimulatorSection simulatorSection_;
   WorldSection worldSection_;
   ParticlesSection throwInParticlesSection_;
   std::vector<SpeciesTypeSection> speciesTypeSections_;
   std::vector<ReactionRuleSection> reactionRuleSections_;
   CopyNumbersSection copyNumbersSection_;
   ParticlePositionSection particlePositionsSection_;

   SectionBase* current_;
   std::array<std::string, 10> parameters_;

   const std::regex regex_pattern_ = std::regex("^\\[(\\w+)\\](.*)");
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const SimulatorSettings& settings)
{
   stream << settings.getSimulatorSection();
   stream << settings.getWorldSection();;
   for (auto& item : settings.getSpeciesTypeSections()) stream << item;
   for (auto& item : settings.getReactionRuleSections()) stream << item;
   stream << settings.getParticlesSection();
   stream << settings.getCopyNumbersSection();
   stream << settings.getParticlePositionsSection();
   stream << std::endl;
   return stream;
}

// --------------------------------------------------------------------------------------------------------------------------------

inline std::istream& operator>>(std::istream& stream, SimulatorSettings& settings)
{
   std::string line;
   while (std::getline(stream, line))
   {
      if (settings.is_section_heading(line))
         settings.create_section_heading(line);
      else
         settings.add_keypair(line);
   }
   return stream;
}

#endif
