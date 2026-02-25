/**
 * @file CLParser.h
 * @brief Minimal command-line argument parsing utilities.
 *
 * This header defines a small hierarchy of parser classes to simplify
 * reading options from @c main (int argc, char* argv[]). It supports:
 * - Templated options bound to user variables (CLOption<T>)
 * - String options (CLStringOption) that may include spaces
 * - Flags (CLFlag) that simply toggle a boolean when present
 */

#ifndef CLPARSER_H
#define CLPARSER_H

#include <iostream>
#include <sstream>
#include <string>

// -----------------------------------------------------------------------------
// Base parser class
// -----------------------------------------------------------------------------

/// \ingroup driver
/**
 * @brief Abstract base for all command-line option parsers.
 *
 * Holds the option string (e.g. "-th" or "--help") and exposes a
 * virtual @c Parse() method to be implemented by derived classes.
 */
class CLParser {
 protected:
  std::string option; //!< Option keyword to match against argv[]

 public:
  /**
   * @brief Construct with an option token.
   * @param opt Option keyword (e.g., "-th", "--help").
   */
  CLParser(std::string opt) : option(opt){};

  /**
   * @brief Construct with an option token.
   * @param opt Option keyword (e.g., "-th", "--help").
   */
  CLParser(const char *opt) : option(opt){};
  
  /// @brief Virtual destructor.
  virtual ~CLParser(){};

  /**
   * @brief Attempt to parse this option from argv[].
   * @param argc Argument count.
   * @param argv Argument vector.
   * @return 1 if found and parsed, -1 if found but invalid, 0 if not found.
   */
  virtual int Parse(int argc, char *argv[]) = 0;
};

// -----------------------------------------------------------------------------
// Generic option class (expects a typed value)
// -----------------------------------------------------------------------------

/// \ingroup driver
/**
 * @brief Typed option expecting a value (e.g., "-th 8", "-l 2").
 * @tparam T Value type to parse (e.g., int, double, std::string).
 */
template <class T>
class CLOption : public CLParser {
 protected:
  T *target;  //!< Pointer to user variable to be set.

 public:
  /**
   * @brief Construct a typed option bound to a user variable.
   * @param opt Option keyword.
   * @param targ Pointer to target variable.
   */
  CLOption(std::string opt, T *targ) : CLParser(opt), target(targ){};

  /**
   * @brief Construct a typed option bound to a user variable.
   * @param opt Option keyword.
   * @param targ Pointer to target variable.
   */
  CLOption(const char *opt, T *targ) : CLParser(opt), target(targ){};

  /// @brief Destructor.
  virtual ~CLOption(){};

  /**
   * @brief Parse the option and extract a typed value from the next token.
   * @param argc Argument count.
   * @param argv Argument vector.
   * @return 1 if found and set, -1 if conversion failed, 0 if not present.
   */
  int Parse(int argc, char *argv[]) {
    std::stringstream strStream(std::stringstream::in);
    for (int i = 1; i < argc; i++) {
      if (std::string(argv[i]) == option) {
        i++;
        strStream.str(argv[i]);
        strStream >> *target;
        if (strStream.fail()) {
          std::cerr << "Invalid argument \"" << argv[i] << "\" for option "
                    << option << "." << std::endl;
          return -1;
        } else
          return 1;
      }
    }
    return 0;
  };
};

// -----------------------------------------------------------------------------
// String option class
// -----------------------------------------------------------------------------

/// \ingroup driver
/**
 * @brief Option that captures a string argument.
 * Example: "-j input.job" stores "input.job".
 */
class CLStringOption : public CLParser {
 protected:
  std::string *target;  //!< Pointer to user string to be set.

 public:
  /**
   * @brief Construct a string option bound to a user string.
   * @param opt Option keyword.
   * @param targ Pointer to target string.
   */
  CLStringOption(std::string opt, std::string *targ)
      : CLParser(opt), target(targ){};
  
  /**
   * @brief Construct a string option bound to a user string.
   * @param opt Option keyword.
   * @param targ Pointer to target string.
   */
  CLStringOption(const char *opt, std::string *targ)
      : CLParser(opt), target(targ){};

  /// @brief Destructor.
  virtual ~CLStringOption(){};

  /**
   * @brief Parse the option and copy the next token as a string.
   * @param argc Argument count.
   * @param argv Argument vector.
   * @return 1 if found and set, 0 if not present.
   */
  int Parse(int argc, char *argv[]) {
    for (int i = 1; i < argc; i++) {
      if (std::string(argv[i]) == option) {
        i++;
        *target = std::string(argv[i]);
        return 1;
      }
    }
    return 0;
  };
};

// -----------------------------------------------------------------------------
// Flag option class
// -----------------------------------------------------------------------------

/// \ingroup driver
/**
 * @brief Boolean flag option with no argument.
 * Example: "--help" sets the bound bool to true.
 */
class CLFlag : public CLParser {
 protected:
  bool *target; //!< Pointer to user bool to be set true if present.

 public:

  /**
   * @brief Construct a flag option bound to a @c bool.
   * @param opt Option keyword.
   * @param targ Pointer to target bool.
   */
  CLFlag(std::string opt, bool *targ) : CLParser(opt), target(targ){};

  /**
   * @brief Construct a flag option bound to a @c bool.
   * @param opt Option keyword.
   * @param targ Pointer to target bool.
   */
  CLFlag(const char *opt, bool *targ) : CLParser(opt), target(targ){};

  /// @brief Destructor.
  virtual ~CLFlag(){};

  /**
   * @brief Parse the option and set the target bool to true if present.
   * @param argc Argument count.
   * @param argv Argument vector.
   * @return 1 if found and set, 0 if not present.
   */
  int Parse(int argc, char *argv[]) {
    for (int i = 1; i < argc; i++) {
      if (std::string(argv[i]) == option) {
        *target = true;
        return 1;
      }
    }
    return 0;
  };
};

#endif
