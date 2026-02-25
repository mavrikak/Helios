/**
 * @file IncidentField.cpp
 * @brief Minimal implementation pieces for @ref IncidentField.
 */

#include "IncidentField.h"

//------------------------------------------------------------------------------
// Lifetime control
//------------------------------------------------------------------------------
IncidentField::~IncidentField() {}

//------------------------------------------------------------------------------
// Global accuracy toggle (static)
//------------------------------------------------------------------------------
bool IncidentField::NeedsAccurate = false;

void IncidentField::EnableAccurate() {
  std::cout << "Using high-accuracy integration for excitation." << std::endl;
  NeedsAccurate = true;
}
