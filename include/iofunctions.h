/**
 * @file iofunctions.h
 * @brief Small helper templates for reading whitespace-separated values from text files
 * while skipping line comments.
 *
 * Provides three simple, header-only utilities:
 * - ReadCommented: reads a token into a typed variable, skipping lines that start with '#'.
 * - ReadCommentedG: same as above but treats '$' as the comment introducer.
 * - ReadAndPrint: convenience for debugging; reads one value and prints it.
 *
 * Each reader consumes input tokens separated by whitespace. If the first non-
 * whitespace character of the next token is the comment character, the rest of
 * that line is skipped. Multiple values per line are supported as long as they
 * are separated by whitespace and appear before any comment token.
 *
 * @note These functions expect a valid, open std::ifstream. They do not alter
 * stream formatting flags.
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// Only include this once during compiling
#ifndef IOFUNCTIONS_H
#define IOFUNCTIONS_H

/// \ingroup driver
/**
 * @brief Read a value of type @p T from a text stream, skipping lines beginning with '#'.
 *
 * Attempts to extract the next non-comment token from @p inStream and parse it
 * into @p target using operator>> on a temporary std::stringstream. If the next
 * token starts with '#', the remainder of that line is discarded and the search
 * continues. The function returns immediately after parsing a single value or
 * if end-of-file/error is reached.
 *
 * @tparam T Destination type (must be stream-extractable via operator>>).
 * @param[in,out] inStream Pointer to an open std::ifstream to read from.
 * @param[out] target Pointer to the destination variable to receive the value.
 * @return 0 on success; 1 on failure to parse or on reaching EOF without reading a value.
 *
 * @par Behavior
 * - Skips tokens that begin with '#', consuming the rest of that line.
 * - Accepts multiple values per line; stops after extracting the first valid value.
 * - Leaves the stream positioned after the consumed token (or end of file on failure).
 *
 * @warning This function assumes that comment markers appear as the first
 * character of a token. Comments that appear *after* a value on the
 * same line must be separated by whitespace (e.g., "42 # note").
 */
template <class T>
int ReadCommented(std::ifstream *inStream, T *target) {
  std::string strRead;
  std::stringstream strStream(std::stringstream::in);

  while ((*inStream)) {
    *inStream >> strRead;
    if (strRead[0] == '#') {
      std::getline(*inStream, strRead);
    } else if ((*inStream)) {
      strStream.str(strRead);
      if (!(strStream >> *target))
        return 1;
      else
        return 0;
    }
  }
  return 1;
}

/// \ingroup driver
/**
 * @brief Read a value of type @p T from a text stream, skipping lines beginning with '$'.
 *
 * Identical semantics to @ref ReadCommented except that lines whose first token
 * starts with '$' are treated as comments.
 *
 * @tparam T Destination type (must be stream-extractable via operator>>).
 * @param[in,out] inStream Pointer to an open std::ifstream to read from.
 * @param[out] target Pointer to the destination variable to receive the value.
 * @return 0 on success; 1 on failure to parse or on reaching EOF without reading a value.
 */
template <class T>
int ReadCommentedG(std::ifstream *inStream, T *target) {
  std::string strRead;
  std::stringstream strStream(std::stringstream::in);

  while ((*inStream)) {
    *inStream >> strRead;
    if (strRead[0] == '$') {
      std::getline(*inStream, strRead);
    } else if ((*inStream)) {
      strStream.str(strRead);
      if (!(strStream >> *target))
        return 1;
      else
        return 0;
    }
  }
  return 1;
}

/// \ingroup driver
/**
 * @brief Debug helper: read a value of type @p T and print it to std::cout.
 *
 * Uses @ref ReadCommented to attempt a read from @p inStream; prints "Error!"
 * if no value could be read/parsed, otherwise prints the value followed by a
 * newline. The function always returns 0 (it is intended for quick debugging).
 *
 * @tparam T Destination type (must be stream-extractable via operator>>).
 * @param[in,out] inStream Pointer to an open std::ifstream to read from.
 * @return Always 0.
 */
template <class T>
int ReadAndPrint(std::ifstream *inStream) {
  T tTmp;
  if (ReadCommented<T>(inStream, &tTmp) != 0)
    std::cout << "Error!" << std::endl;
  else
    std::cout << tTmp << std::endl;
  return 0;
}

#endif

