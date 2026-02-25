/**
 * @file SIEForm.cpp
 * @brief Implementation for the SIEForm class.
 */

#include "SIEForm.h"

// Construct SIEForm and compute angular frequency from vacuum wavelength.
SIEForm::SIEForm(std::vector<RWGFun> *inRWGFuns, GreenF *inGrnFun,
                 std::vector<IncidentField *> *inIncFields,
                 dcmplx inVacWavelength)
    : rwgFuns(inRWGFuns),
      grnFun(inGrnFun),
      incFields(inIncFields),
      omega(2 * PI / inVacWavelength * CVAC) {}

// Virtual destructor (no ownership implied).
SIEForm::~SIEForm() {}
