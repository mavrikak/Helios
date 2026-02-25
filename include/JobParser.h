/**
 * @file JobParser.h
 * @brief Lightweight XML-like tag parser for the solver job file.
 *
 * Provides a minimal, streaming parser that understands a restricted, whitespace-
 * separated tag syntax, e.g.:
 *
 * @f$<job> <mesh> mesh.mesh </mesh> <excitation> ... </excitation> </job>@f$
 *
 * The class keeps a stack of open tags to allow nested blocks and exposes
 * helpers to read typed values and to jump to the end of the current tag.
 */

#ifndef JOBPARSER_H
#define JOBPARSER_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

/// \ingroup driver
/**
 * @class JobParser
 * @brief Streaming, stack-based reader for simple tag-delimited job files.
 *
 * The parser operates directly on an std::ifstream. It tracks the current open
 * tag using an internal stack ( @p tags ). Opening tags push onto the stack;
 * closing tags pop it. The API exposes methods to open/close a file, read the
 * next tag token, query the current tag, read a typed value, and skip to the end
 * of the current tag block.
 *
 * @warning This is **not** a full XML parser. It expects tokens separated by
 * whitespace, with tags like `<name>` and `</name>` appearing as single
 * tokens. Attributes are not supported.
 */
class JobParser {
 private:
  std::ifstream reader;           //!< Underlying input stream
  std::vector<std::string> tags;  //!< Stack of open tags
  std::string buffer;             //!< Scratch buffer
  int nextChar();                 //!< Skip whitespace (ws); position on next non-ws
  std::string tag;                //!< Most recently read tag token

 public:
  /** @name Construction & file control */
  //{@
  /** 
   * @brief Default-construct an empty parser (no file opened). 
   */
  JobParser();

  /** 
   * @brief Construct and open a job file; prints an error on failure. 
   * @param jobFile The name of the job file to open.
   */
  JobParser(std::string jobFile);

  /** @brief Open a job file.
   *  @param jobFile The name of the job file to open.
   *  @return 0 on success, 1 on failure.
   */
  int open(std::string jobFile);

  /** 
   * @brief Close the currently open file. Always returns 0.
   * @return Always returns 0.
   */
  int close();
  //}@

  /** @name Tag navigation */
  //{@
  /**
  * @brief Read the next tag token and update the open-tag stack.
  *
  * If the next non-whitespace character is '<', the token following it is
  * interpreted as a tag (with the trailing '>' stripped). Opening tags are
  * pushed to the stack; closing tags (starting with '/') pop if they match the
  * current top.
  *
  * @return The name of the current open tag after the update, or an empty
  * string if there are no open tags.
  */
  std::string readTag();

  /** 
   * @brief Return the current open tag (stack top), or empty string if none. 
   * @return The name of the current open tag after the update, or an empty
   * string if there are no open tags.
   */
  std::string getTag();

  /** 
   * @brief Return the most recently read raw tag token (without '>'). 
   * @return The most recently read tag token (raw, without trailing '>').
   */
  std::string lastRead();

  /**
  * @brief Skip tokens until the matching closing tag of the current block.
  *
  * Reads tokens using the stream extraction operator until `</getTag()>` is
  * encountered, then returns. Intended to quickly skip over blocks you do not
  * wish to parse.
  * @return 0 on completion.
  */
  int endOfTag();


  /** 
   * @brief Seek the stream back to the beginning. 
   * @return 0 on success.
   */
  int moveTop();
  //}@
  
  /** @name Value extraction */
  //{@
  /** 
   * @brief True if there are no open tags on the stack. 
   * @return True if there are no open tags on the stack.
   */
  bool empty();

  /**
  * @brief Read a value of type @p Type from the stream using operator>>.
  * @tparam Type to read.
  * @return The value read.
  *
  * Consumes the next whitespace-delimited token from the underlying stream and
  * converts it to @p Type. No comment handling or tag awareness is performed
  * here; prefer to call this after @ref readTag has positioned you inside the
  * right block.
  */
  template <class Type>
  Type read();
  //}@
};

// ------------------- Inline template implementation -------------------------
/**
 * @brief Read a value of type @p Type from the stream using operator>>.
 * @tparam Type to read.
 * @return The value read.
 *
 * Consumes the next whitespace-delimited token from the underlying stream and
 * converts it to @p Type. No comment handling or tag awareness is performed
 * here; prefer to call this after @ref readTag has positioned you inside the
 * right block.
 */
template <class Type>
Type JobParser::read() {
  Type t;
  reader >> t;  // Assumes Type supports operator>>(istream&, Type&)
  return t;
}
#endif
