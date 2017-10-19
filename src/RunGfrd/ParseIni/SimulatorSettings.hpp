#ifndef SIMULATORSETTINGS_HPP
#define SIMULATORSETTINGS_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <vector>
#include "ParserExceptions.hpp"
#include "VariablesSection.hpp"
#include "SimulatorSection.hpp"
#include "WorldSection.hpp"
#include "SpeciesTypeSection.hpp"
#include "ReactionRuleSection.hpp"
#include "ParticlesSection.hpp"
#include "CopyNumbersSection.hpp"
#include "ParticlePositionsSection.hpp"
#include "ReactionRecordSection.hpp"
#include "ProgressSection.hpp"

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
      bool is_match = std::regex_search(line, match, regex_section_);
      THROW_UNLESS_MSG(illegal_section, is_match, line);

      auto heading = match[1].str();
      create_section(heading);
      current_->set_vars(&variablesSection_);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool is_section_heading(const std::string& line) const
   {
      return std::regex_match(line, regex_section_);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void add_keypair(const std::string& line) const
   {
      if (line.empty()) return;
      
      auto key = std::string();
      auto value = std::string();
      bool is_match = get_name_value(line, key, value);
      if (is_match)
      {
         THROW_UNLESS_MSG(illegal_section_key, current_, "Key without section: " << line)
         current_->set_keypair(key, value);
      }
      else
      {
         std::smatch match;
         auto regex = std::regex("^\\s*([\\w<>\\$\\(\\)\\^\\']+[^;]*).*");
         bool is_match = std::regex_match(line, match, regex);
         THROW_UNLESS_MSG(illegal_section_value, !is_match, "Invalid line encountered: " << line)
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void add_variable(const std::string& line)
   {
      auto key = std::string();
      auto value = std::string();
      bool is_match = get_name_value(line, key, value);
      if (is_match) variablesSection_.add_variable(key, value, true);
      else THROW_EXCEPTION(illegal_argument, "Unexpected assignment in '" << line << "'.")
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:

   bool get_name_value(const std::string& line, std::string& name, std::string& value) const
   {
      std::smatch match;
      auto regex = std::regex("^\\s*(\\w+)\\s*=\\s*([\\w<>\\$\\(\\)\\^\\']+[^;]*).*");
      bool is_match = std::regex_search(line, match, regex);
      name = match[1].str();
      value = match[2].str();

      auto iend = value.rbegin();
      size_t c = 0;
      while (iend != value.rend() && (*iend == ' ' || *iend == '\t'))  ++c,++iend;

      value = value.substr(0, value.length() - c);
      return is_match;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void create_section(const std::string& heading)
   {
      if (heading == VariablesSection::section_name()) { current_ = &variablesSection_; return; }
      if (heading == SimulatorSection::section_name()) { current_ = &simulatorSection_; return; }
      if (heading == WorldSection::section_name()) { current_ = &worldSection_; return; }
      if (heading == ParticlesSection::section_name()) { particlesSections_.emplace_back(); current_ = &(*particlesSections_.rbegin()); return; }
      
      if (heading == CopyNumbersSection::section_name()) { if (copyNumbersSection_.get() == nullptr) copyNumbersSection_ = std::make_unique<CopyNumbersSection>(); current_ = copyNumbersSection_.get(); return; }
      if (heading == ParticlePositionSection::section_name()) { if (particlePositionsSection_.get() == nullptr) particlePositionsSection_ = std::make_unique<ParticlePositionSection>(); current_ = particlePositionsSection_.get(); return; }
      if (heading == ReactionRecordSection::section_name()) { if (reactionRecordSection_.get() == nullptr) reactionRecordSection_ = std::make_unique<ReactionRecordSection>(); current_ = reactionRecordSection_.get(); return; }
      if (heading == ProgressSection::section_name()) { if (progressSection_.get() == nullptr) progressSection_ = std::make_unique<ProgressSection>(); current_ = progressSection_.get(); return; }
      
      if (heading == ReactionRuleSection::section_name()) { reactionRuleSections_.emplace_back(); current_ = &(*reactionRuleSections_.rbegin()); return; }
      if (heading == SpeciesTypeSection::section_name()) { speciesTypeSections_.emplace_back(); current_ = &(*speciesTypeSections_.rbegin()); return; }
      
      THROW_EXCEPTION(illegal_section, "Section '" << heading << "' not recognized!");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

public:

   VariablesSection& getVariablesSection() { return variablesSection_; }
   const VariablesSection& getVariablesSection() const { return variablesSection_; }
   const SimulatorSection& getSimulatorSection() const { return simulatorSection_; }
   const WorldSection& getWorldSection() const { return worldSection_; }
   const std::vector<SpeciesTypeSection>& getSpeciesTypeSections() const { return speciesTypeSections_; }
   const std::vector<ReactionRuleSection>& getReactionRuleSections() const { return reactionRuleSections_; }
   const std::vector<ParticlesSection>& getParticlesSections() const { return particlesSections_; }
   const CopyNumbersSection* getCopyNumbersSection() const { return copyNumbersSection_.get(); }
   const ParticlePositionSection* getParticlePositionsSection() const { return particlePositionsSection_.get(); }
   const ReactionRecordSection* getReactionRecordSection() const { return reactionRecordSection_.get(); }
   const ProgressSection* getProgressSection() const { return progressSection_.get(); }

   // --------------------------------------------------------------------------------------------------------------------------------

private:

   VariablesSection variablesSection_;
   SimulatorSection simulatorSection_;
   WorldSection worldSection_;
   std::vector<ParticlesSection> particlesSections_;
   std::vector<SpeciesTypeSection> speciesTypeSections_;
   std::vector<ReactionRuleSection> reactionRuleSections_;
   std::unique_ptr<CopyNumbersSection> copyNumbersSection_;
   std::unique_ptr<ParticlePositionSection> particlePositionsSection_;
   std::unique_ptr<ReactionRecordSection> reactionRecordSection_;
   std::unique_ptr<ProgressSection> progressSection_;

   SectionBase* current_;
   const std::regex regex_section_ = std::regex("^\\s*\\[\\s*(\\w+)\\s*\\](.*)");
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const SimulatorSettings& settings)
{
   stream << settings.getVariablesSection();
   stream << settings.getSimulatorSection();
   stream << settings.getWorldSection();
   for (auto& item : settings.getSpeciesTypeSections()) stream << item;
   for (auto& item : settings.getReactionRuleSections()) stream << item;
   for (auto& item : settings.getParticlesSections()) stream << item;
   stream << settings.getCopyNumbersSection();
   stream << settings.getParticlePositionsSection();
   stream << settings.getReactionRecordSection();
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
