#include "pamCPI.h"

using namespace std;

pamCPI::pamCPI(prefAttachModel *pa, const int projStep, const int collectInterval, const int offManifoldWait, const int nMicroSteps) : projStep(projStep), collectInterval(collectInterval), offManifoldWait(offManifoldWait), nMicroSteps(nMicroSteps) {
    pam = pa;
};

void pamCPI::runCPI(const int nSteps) {
    //after waiting for the system to reach the slow manifold, collect data every collectInterval number of steps
    int totalSteps = 0;
    int saveDataInterval = nMicroSteps/100;
    while(totalSteps < nSteps) {
	int microStep;
	vector<graphData> toPlot;
	vector<graphData> toProject;
	for(microStep = 0; microStep < nMicroSteps; microStep++) {
	    if(microStep < offManifoldWait) {
		if((microStep+1)%saveDataInterval == 0) {
		    toPlot.push_back(pam->step(true));
		}
		else {
		    pam->step(false);
		}
	    }
	    else {
		int nOnManifoldSteps = microStep - offManifoldWait;
		if((nOnManifoldSteps+1)%collectInterval == 0) {
		    if((microStep+1)%saveDataInterval == 0) {
			graphData d = pam->step(true);
			toPlot.push_back(d);
			toProject.push_back(d);
		    }
		    else {
			toProject.push_back(pam->step(true));
		    }
		}
		else {
		    if((microStep+1)%saveDataInterval == 0) {
			toPlot.push_back(pam->step(true));
		    }
		    else {
			toProject.push_back(pam->step(false));
		    }
		}
	    }
	}
	//project collected data forward, save toPlot data, clear all vectors, update totalSteps
    }
}
		    
