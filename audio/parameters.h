#ifndef EEGFILTER_PARAMETERS_H
#define EEGFILTER_PARAMETERS_H

#include "dnf.h"

const char* const p300Path = "../../gla_researchdata_1258/EEG_recordings/participant%03d/rawp300.tsv";
const char* const tasksPath = "../../gla_researchdata_1258/EEG_recordings/participant%03d/%s.tsv";

const Neuron::actMethod ACTIVATION = Neuron::Act_Tanh; //keep default activation method

// number of subjects
const int nSubj = 20; // audio in subj no.5

// pre-filtering
const int filterorder = 2;
const double innerHighpassCutOff = 20000; // Hz 0.5 400
const double outerHighpassCutOff = 400; // Hz 5, 20000
const double LaplaceCutOff = 20000; // Hz 0.5 , 20000

const double powerlineFrequ = 20000; // Hz 50 , 20000
const double bsBandwidth = 10000; // Hz 2.5, 10000

//creat circular buffers for plotting
const int bufferLength = 1000 ;

// dnf learning rate
const double dnf_learning_rate_p300 = 10; //10
const double dnf_learning_rate_tasks = 2.5; //2.5

// dnf number of layers
const int NLAYERS = 6;

// LMS learning rates
const double lms_learning_rate_p300 = 0.05; //0.04, was 0.005 in lms data
const double lms_learning_rate_tasks = 0.01; //0.01

const double inner_gain = 10; //changed from 1000 to 1, SN
const double outer_gain = 1; //changed from 1000 to 1, N
const double remover_gain = 1;

// Very slow
// #define SAVE_WEIGHTS

#endif //EEGFILTER_PARAMETERS_H
