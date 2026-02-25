/**
 * @file JobParser.cpp
 * @brief Implementation of the lightweight tag parser used for job files.
 */

#include "JobParser.h"
#include "iofunctions.h"

//------------------------------------------------------------------------------
// Consrtuctors
//------------------------------------------------------------------------------
/** @brief Default constructor: no file opened. */
JobParser::JobParser() {}

/**
 * @brief Construct and attempt to open @p jobFile.
 * @details Prints an error message if the file cannot be opened; leaves the
 * stream closed in that case.
 */
JobParser::JobParser(std::string jobFile) {
  reader.open(jobFile.c_str());
  if (!reader.is_open()) {
    std::cout << "Error opening job file: " << jobFile << std::endl;
    reader.close();
  }
}

//------------------------------------------------------------------------------
// File control
//------------------------------------------------------------------------------
/** @brief Close the underlying stream. */
int JobParser::close() {
  reader.close();
  return 0;
}

/** @brief True if there are no open tags. */
bool JobParser::empty() {
  if (tags.size() == 0) {
    return true;
  } else {
    return false;
  }
}

/**
 * @brief Open @p jobFile.
 * @return 0 on success; 1 on failure (and prints an error message).
 */
int JobParser::open(std::string jobFile) {
  reader.open(jobFile.c_str());
  if (!reader.is_open()) {
    std::cout << "Error opening job file: " << jobFile << std::endl;
    reader.close();
    return 1;
  } else {
    return 0;
  }
}

//------------------------------------------------------------------------------
// Tag navigation
//------------------------------------------------------------------------------
/**
 * @brief Read the next tag (opening or closing) and update the tag stack.
 *
 * Inspects the next non-whitespace character; if it is '<', extracts the tag
 * token (stripping the trailing '>'). If the tag does not start with '/', it is
 * pushed onto the tag stack. Otherwise, the stack is popped if the closing tag
 * matches the current top; if it does not, a warning is printed.
 *
 * @return The name of the currently open tag (stack top) after processing, or
 * an empty string if the stack becomes empty.
 */
std::string JobParser::readTag() {
  char id;
  nextChar();
  reader.get(id);
  reader.unget();
  if (id == '<') {
    reader.get(id);
    reader >> tag;
    tag = tag.replace(tag.size() - 1, 1, "");
    if (tag[0] != '/') {
      tags.push_back(tag);
    } else {
      std::string comp(tag);
      if (comp.replace(0, 1, "") == JobParser::getTag()) {
        tags.pop_back();
      } else {
        std::cout << "Incorrect use of tags in job file" << std::endl;
      }
    }
  }
  if (tags.size() != 0) {
    return tags.back();
  }

  return "";
}

/** @brief The most recently read tag token (raw, without trailing '>'). */
std::string JobParser::lastRead() { return tag; }


/**
 * @brief Advance the stream to the next non-whitespace character and back up one.
 * @return Always 0.
 */
int JobParser::nextChar() {
  char id;
  reader.get(id);
  while (((id == ' ') || (id == '\n') || (id == '\r') || (id == '\t')) &&
         (!reader.eof())) {
    reader.get(id);
  }
  reader.unget();
  return 0;
}

/** @brief Return the current open tag (stack top), or empty string if none. */
std::string JobParser::getTag() {
  if (tags.size() != 0) {
    return tags.back();
  } else {
    return "";
  }
}

//------------------------------------------------------------------------------
// Block helpers
//------------------------------------------------------------------------------


/**
 * @brief Consume tokens until the closing tag of the current block is found.
 * @details Builds the token `</getTag()>` and keeps reading whitespace-delimited
 * tokens until an exact match is encountered.
 * @return 0 when the closing tag token has been read (or EOF reached).
 */
int JobParser::endOfTag() {
  std::string endTag("</" + getTag() + ">");
  std::string word;
  while (reader >> word && word != endTag) {
  }
  return 0;
}

/** @brief Seek to the beginning of the file (sets get pointer to beginning). */
int JobParser::moveTop() {
  reader.seekg(0, reader.beg);
  return 0;
}
