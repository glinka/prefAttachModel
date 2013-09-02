#include <algorithm>
#include "pamCPI.h"
#include "calcGraphProps.h"
#include "fitCurves.h"

using namespace std;

pamCPI::pamCPI(const int n, const int m, const double kappa, const int projStep, const int collectInterval, const int offManifoldWait, const int nMicroSteps) : prefAttachModel(n, m, kappa), projStep(projStep), collectInterval(collectInterval), offManifoldWait(offManifoldWait), nMicroSteps(nMicroSteps) {
};

void pamCPI::runCPI(const int nSteps) {
    initGraph();
    //set up save file and write header
    vector<double> forFile;
    forFile.push_back(projStep);
    forFile.push_back(offManifoldWait);
    forFile.push_back(nMicroSteps);
    vector<string> forFileStrs;
    forFileStrs.push_back("projStep");
    forFileStrs.push_back("offManifoldWait");
    forFileStrs.push_back("nMicroSteps");
    ofstream *paDataCPI = createFile("paDataCPI", forFile, forFileStrs);
    ofstream *projData = createFile("projData", forFile, forFileStrs);
    //after waiting for the system to reach the slow manifold, collect data every collectInterval number of steps
    int totalSteps = 0;
    int saveDataInterval = nMicroSteps/100;
    while(totalSteps < nSteps) {
	int microStep;
	vector<graphData> toPlot;
	vector<graphData> toProject;
	vector<double> time;
	for(microStep = 0; microStep < nMicroSteps; microStep++) {
	    if(microStep < offManifoldWait) {
		if((microStep+1)%saveDataInterval == 0) {
		    toPlot.push_back(step(true));
		}
		else {
		    step(false);
		}
	    }
	    else {
		int nOnManifoldSteps = microStep - offManifoldWait;
		if((nOnManifoldSteps+1)%collectInterval == 0) {
		    time.push_back(totalSteps);
		    if((microStep+1)%saveDataInterval == 0) {
			graphData d = step(true);
			toPlot.push_back(d);
			toProject.push_back(d);
		    }
		    else {
			toProject.push_back(step(true));
		    }
		}
		else {
		    if((microStep+1)%saveDataInterval == 0) {
			toPlot.push_back(step(true));
		    }
		    else {
			step(false);
		    }
		}
	    }
	    totalSteps++;
	}
	//project collected data forward, save toPlot data, clear all vectors, update totalSteps
	vector<vector<vector<double> > > data;
	vector<vector<double> > mat;
	vector<double> v;
	for(vector<graphData>::iterator d = toProject.begin(); d != toProject.end(); d++) {
	    data.push_back(mat);
	    int i, j;
	    for(i = 0; i < n; i++) {
		data.back().push_back(v);
		for(j = 0; j < n; j++) {
		    (data.back()[i]).push_back((*d).A[i][j]);
		}
	    }
	}
	project(data, time, *projData);
	int nPlotPts = toPlot.size();
	//hooray for vlas
	graphData toPlotAry[nPlotPts];
	int i;
	for(i = 0; i < nPlotPts; i++) {
	    toPlotAry[i] = toPlot[i];
	}
	saveData(toPlotAry, nPlotPts, *paDataCPI);
	toPlot.clear();
	toProject.clear();
	totalSteps += projStep;
    }
    paDataCPI->close();
    delete paDataCPI;
}
/**
******************** TODO ********************
takes adjacency matrices vs time in the form of a vector of matrices (vector of a vector of a vector) and simple vector, projects with some funciton
would be nice to be able to specify the function, but would probably require rewriting many of calcGraphProps.cc code to take vectors, or, better, to 
define throughout:
typedef vector<vector<double>> graph;
and pass all calcGraphProps fns a 'graph'
**********************************************
**/
void pamCPI::project(vector<vector<vector<double> > > &data, vector<double> &time, ofstream &projData) {
    //assert data.size() == time.size()
    int nPts = time.size();
    int i, j, k;
    vector<vector<double> > eigVectFittedCoeffs;
    vector<double> line;
    //fit eigvect with fifth order polynomial
    fxs toFitEigVects;
    fx constOffset = [] (double x) -> double { 
	return 1;
    };
    toFitEigVects.push_back(constOffset);
    toFitEigVects.push_back([] (double x) { return x;});
    toFitEigVects.push_back([] (double x) { return x*x;});
    toFitEigVects.push_back([] (double x) { return x*x*x;});
    toFitEigVects.push_back([] (double x) { return x*x*x*x;});
    toFitEigVects.push_back([] (double x) { return x*x*x*x*x;});
    //fit coeffs evolution over time with line
    fxs toFitCoeffs;
    toFitCoeffs.push_back(constOffset);
    toFitCoeffs.push_back([] (double x) { return x;});
    for(i = 0; i < n; i++) {
	line.push_back(i);
    }
    for(i = 0; i < nPts; i++) {
	vector<vector<double> > currentAdjMat = data[i];
	//want last column, corresponds to eigenvector with largest eigenvalue
	int **tempA = new int*[n];
	for(j = 0; j < n; j++) {
	    tempA[j] = new int[n];
	    for(k = 0; k < n; k++) {
		tempA[j][k] = currentAdjMat[j][k];
	    }
	}
	double *eigVals = calcGraphProps::getAdjEigVals(tempA, n);
	double **eigVects = calcGraphProps::getAdjEigVects(tempA, n);
	double maxEigVal = eigVals[n-1];
	vector<double> leadingEigVect;
	for(j = 0; j < n; j++) {
	    leadingEigVect.push_back(eigVects[j][n-1]);
	}
	for(j = 0; j < n; j++) {
	  delete[] eigVects[j];
	}
	sort(leadingEigVect.begin(), leadingEigVect.end());
	delete eigVects;
	delete[] eigVals;
	eigVectFittedCoeffs.push_back(fitCurves::fitFx(line, leadingEigVect, toFitEigVects));
	eigVectFittedCoeffs.back().push_back(maxEigVal);
	for(j = 0; j < n; j++) {
	    delete[] tempA[j];
	}
	delete tempA;
    }
    int nCoeffs = eigVectFittedCoeffs[0].size();
    vector<vector<double> > reshapedCoeffs;
    vector<double> v;
    //coeffs used to fit time evolution of (coefficients used to fit time evolution of eigenvectors)
    vector<vector<double> > timeEvoCoeffs;
    for(i = 0; i < nCoeffs; i++) {
	reshapedCoeffs.push_back(v);
	for(j = 0; j < nPts; j++) {
	    reshapedCoeffs[i].push_back(eigVectFittedCoeffs[j][i]);
	}
	timeEvoCoeffs.push_back(fitCurves::fitFx(time, reshapedCoeffs[i], toFitCoeffs));
    }
    int newTime = time.back() + projStep;
    int nSecondaryCoeffs = toFitCoeffs.size();
    vector<double> newCoeffs;
    /**************************************** TODO ****************************************
decrease time vector to be the same during each projection, else values will become overly large
    **************************************** TODO ****************************************/
    for(i = 0; i < nCoeffs; i++) {
	double eval = 0;
	for(j = 0; j < nSecondaryCoeffs; j++) {
	    eval += (*toFitCoeffs[j])(newTime)*timeEvoCoeffs[i][j];
	}
	newCoeffs.push_back(eval);
    }
    double newEigVal = newCoeffs.back();
    newCoeffs.pop_back();
    nCoeffs = newCoeffs.size();
    vector<double> newEigVect;
    for(i = 0; i < n; i++) {
	double eval = 0;
	//eigval appended to end of eigVectFittedCoeffs, iterate only to nCoeffs - 1 as last newCoeff is eigVal
	for(j = 0; j < nCoeffs; j++) {
	    eval += (*toFitEigVects[j])(line[i])*newCoeffs[j];
	}
	newEigVect.push_back(eval);
    }
    int **newA = new int*[n];
    vector<vector<double> > toSaveRecon;
    for(i = 0; i < n; i++) {
	newA[i] = new int[n];
	toSaveRecon.push_back(v);
	for(j = 0; j < n; j++) {
	    //need to round to nearest even on diag
	    if(i == j) {
		int remainder = ((int) newEigVal*newEigVect[i]*newEigVect[j]);
		remainder = remainder % 2;
		newA[i][j] = ((int) newEigVal*newEigVect[i]*newEigVect[j]) + remainder;
	    }
	    else {
		newA[i][j] = (int) (newEigVal*newEigVect[i]*newEigVect[j] + 0.5);
	    }
	    toSaveRecon[i].push_back(newA[i][j]);
	}
    }
    initGraph(newA);
    for(i = 0; i < n; i++) {
	delete[] newA[i];
    }
    delete newA;
    //save data for visual comparison of reconstructions:
    vector<vector<double> > toSavePreRecon = data.back();
    saveData(toSavePreRecon, projData);
    saveData(toSaveRecon, projData);
}
