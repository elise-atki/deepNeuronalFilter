// Usage Examples
//
// LMS demo
// removes 50Hz from an Audio signal with a powerline reference
// here generated with a simple sine. Should be ideally
// also measured with another ADC channel and fed into
// the filter.

// This is the only include you need
#include <Fir1.h>
#include <iostream>
#include <iomanip>

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>

#define NTAPS 500
#define LEARNING_RATE 0.005

int main (int,char**)
{
	// inits the filter
	Fir1 fir(NTAPS);
	fir.setLearningRate(LEARNING_RATE);

	FILE *finput = fopen("recording2voxs+n.dat","rt");
	FILE *foutput = fopen("recfilteredout.dat","wt");
	FILE *finputref = fopen("recording2bassn.dat","rt");
      
	
	for(int i=0;i<1510000;i++) //i<no. of samples in voxs+n
	{
	  float input_signal, timestamp1;		
	  fscanf(finput,"%f %f\n",&timestamp1, &input_signal);
	  
	  float ref_noise, timestamp2;
	  fscanf(finputref, "%f %f/n",&timestamp2, &ref_noise);  
	  //float ref_noise = sin(2*M_PI/20*i);
	  
	  float canceller = fir.filter(ref_noise);
	  float output_signal = input_signal - canceller;
    

	  
	  fir.lms_update(output_signal); 
	  fprintf(foutput,"%f %f %f\n",output_signal,canceller,ref_noise);
	}
	
	fclose(finput);
	fclose(foutput);
	fclose(finputref);
	fprintf(stderr,"Written the filtered ECG to 'recfilteredout.dat'\n");
}



/////////////////////////


// This is the only include you need
#include <Fir1.h>
#include <iostream>
#include <iomanip>
#include "dnf.h"
#include "parameters.h"
#include <fstream>
#include <string>
#include <stdio.h>
#include <boost/circular_buffer.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/lexical_cast.hpp>
#include <Iir.h>
#include <chrono>
#include <string>
#include <ctime>
#include <thread>         // std::thread
#include <memory>
#include <numeric>

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>

// take from http://stackoverflow.com/a/236803/248823
void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

int main(int,char**,const char** taskinit=nullptr)
{

  std::string outpPrefix = "audioout";
  //initalising constants
  int fs = 44100;
  if (nullptr!=taskinit){
    fs = 44100;
    outpPrefix = taskinit;
  }
  
  const int samplesNoLearning = 3 * fs / innerHighpassCutOff;//iphcofrom parameters.h
  const int nTapsDNF = 100;

  //do i need to add in buffers here?

  fstream dnf_file; //opens up data files
  fstream inner_file;
  fstream outer_file;
  fstream lms_file;
  fstream laplace_file;
  ifstream data_infile; //reads in data
  fstream wdistance_file;
  //! FILE * data_infile = fopen("rawp300.tsv", "rt"); // data_infile.open(file);
#ifdef SAVE_WEIGHTS
  fstream weight_file;
#endif


  long count = 0;

  //no interactive window here
  
  //begin DNF

  DNF dnf(NLAYERS,nTapsDNF,fs,ACTIVATION); 

  boost::circular_buffer<double> intrig_delay(dnf.getSignalDelaySteps());// where does getsigdelay come from - error!?

  //creates outptfiles
  dnf_file.open(outpPrefix + "/dnf.tsv", fstream::out);//outpPrefix+"/subject" + sbjct + 
  inner_file.open(outpPrefix +"/inner.tsv", fstream::out);
  outer_file.open(outpPrefix +"/outer.tsv", fstream::out);
  lms_file.open(outpPrefix +"/lms.tsv", fstream::out);
  laplace_file.open(outpPrefix +"/laplace.tsv", fstream::out);

#ifdef SAVE_WEIGHTS
  weight_file.open(outpPrefix+"/lWeights.tsv", fstream::out);
#endif
  wdistance_file.open(outpPrefix+"/weight_distance.tsv", fstream::out);

  //what is fullpath2data for in other code?
  
  //do i need the filters here?

  Fir1 lms_filter(nTapsDNF);

  while (!data_infile.eof()){ //for i<< small no. samples to test
    count++;

    //initialising signals
    
    double podcastsn = 0;
    double bassn = 0;
    double p300trig = 0; //what is function of this? (used to pick subject?)


    if (nullptr == taskinit {
	data_infile >> bassn >> podcastsn >> p300trig; //opens data_infile
      } else {
	    
    
	      std::string line;
	      std::getline(data_infile,line);

	      vector<string> row_values;
	      split(line, '\t', row_values); //ERROR with split, why? - void1 solved?
	      if (row_values.size()>7) {
		podcastsn = boost::lexical_cast<double>(row_values[7]); //converts string input to double no. val, why 7&8 (1&2?)
		bassn = boost::lexical_cast<double>(row_values[8]);
	      }
	      p300trig = 0;
      }
      

  //adjustment & amplification (leave for now, include HPF later) filtering
  
  //weight/learning rate setup
    double f_nn = dnf.filter(podcastsn,bassn); //podcastsn isn't filtered

    double w_eta = dnf.learningrate_p300;
      if (nullptr != taskinit); {
	w_eta = dnf_learning_rate_tasks;
      }
      if (count > (samplesNoLearning+nTapsDNF)){
	dnf.getNet().setLearningRate(w_eta,0);
      } else {
       	dnf.getNet().setLearningRate(0, 0);
      }
		

  //saving weights to file
#ifdef SAVE_WEIGHTS
      NNO.snapWeights(outpPrefix, "audioout"); 
#endif

      wdistance_file << dnf.getNet().getWeightDistance();
      for(int i=0; i < NLAYERS; i++ ) {
	wdistance_file << "\t" << dnf.getNet().getLayerWeightDistance(i);
      }
      wdistance_file << endl;


      // Do Laplace filter (no filters input yet)
      //double laplace = laplaceHP.filter(podcastsn - bassn); //no filtering
  
      //laplace = laplaceBS.filter(laplace);

      // Do LMS filter (may use my own here?)
      if (count > (samplesNoLearning+nTapsDNF)){
	if (nullptr == taskinit); {
	  lms_filter.setLearningRate(lms_learning_rate_p300);//from parameters.h
	} else {
	  lms_filter.setLearningRate(lms_learning_rate_tasks);
	}
      } else {
	lms_filter.setLearningRate(0);
      }
  
      double corrLMS = lms_filter.filter(podcastsn);
      double lms_output = dnf.getDelayedSignal() - corrLMS;
      if (count > (samplesNoLearning+nTapsDNF)){
	lms_filter.lms_update(lms_output);
      }

      //save signals in files
  
      // undo the gain so that the signal is again in volt, no gain rn so not needed (/inner_gain or /outer_gain in parameters.h)
      //laplace_file << laplace << "\t" << p300trig << endl; 
      inner_file << dnf.getDelayedSignal() << "\t" << endl; //delayed signal through dnf.h ? *
      outer_file << podcastsn << "\t" << endl; //*removed <<delayedp300rtigger as no delay? CAUSING ERROR
      dnf_file << dnf.getOutput() << "\t" << dnf.getRemover() << "\t" << endl;//*
      lms_file << lms_output << "\t" << corrLMS << "\t" << endl;


  dnf.getNet().snapWeights(outpPrefix, "audioout");// what's this line do i need it?
  data_infile.close();
  dnf_file.close();
  inner_file.close();
  outer_file.close();
  lms_file.close();
  laplace_file.close();
  wdistance_file.close();
#ifdef SAVE_WEIGHTS
  weight_file.close();
#endif
}
  /////

  int main (int,char**)
  {
    // inits the filter
    Fir1 fir(NTAPS);
    fir.setLearningRate(LEARNING_RATE);

    FILE *finput = fopen("recording2voxs+n.dat","rt");
    FILE *foutput = fopen("recfilteredout.dat","wt");
    FILE *finputref = fopen("recording2bassn.dat","rt");
      
	
    for(int i=0;i<1510000;i++) //i<no. of samples in voxs+n
      {
	float input_signal, timestamp1;		
	fscanf(finput,"%f %f\n",&timestamp1, &input_signal);
	  
	float ref_noise, timestamp2;
	fscanf(finputref, "%f %f/n",&timestamp2, &ref_noise);  
	//float ref_noise = sin(2*M_PI/20*i);
	  
	float canceller = fir.filter(ref_noise);
	float output_signal = input_signal - canceller;
    

	  
	fir.lms_update(output_signal); 
	fprintf(foutput,"%f %f %f\n",output_signal,canceller,ref_noise);
      }
	
    fclose(finput);
    fclose(foutput);
    fclose(finputref);
    fprintf(stderr,"Written the filtered ECG to 'recfilteredout.dat'\n");
  }
