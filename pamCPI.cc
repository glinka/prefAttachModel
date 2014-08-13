// #include <sys/types.h>
// #include <sys/stat.h>
// #include <unistd.h>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <mpi.h>
#include "pamCPI.h"
#include "calcGraphProps.h"
#include "fitCurves.h"
#include "util_fns.h"

using namespace std;

pamCPI::pamCPI(const int n, const int m, const double kappa, const int projStep, const int collectInterval, const int offManifoldWait, const int nMicroSteps, const int save_interval) : prefAttachModel(n, m, kappa), projStep(projStep), collectInterval(collectInterval), offManifoldWait(offManifoldWait), nMicroSteps(nMicroSteps), save_interval(save_interval) {
};

void pamCPI::runCPI(const int nSteps, const string init_type, const string dir, const string run_id) {
  // TESTING
  // open files once to clear them, afterwards data will be appended
  ofstream times_out("./deg_cpi_data/times" + run_id + ".csv");
  ofstream pre_proj_degs_out("./deg_cpi_data/pre_proj_degs" + run_id + ".csv");
  ofstream post_proj_degs_out("./deg_cpi_data/post_proj_degs" + run_id + ".csv");
  ofstream fitted_coeffs_out("./deg_cpi_data/fitted_coeffs" + run_id + ".csv");
  times_out << "n=" << n << ",proj_step=" << projStep << ",off_manifold_wait=" << offManifoldWait << ",collection_interval=" << collectInterval << ",nms=" << nMicroSteps << endl;
  pre_proj_degs_out << "n=" << n << ",proj_step=" << projStep << ",off_manifold_wait=" << offManifoldWait << ",collection_interval=" << collectInterval << ",nms=" << nMicroSteps << endl;
  post_proj_degs_out << "n=" << n << ",proj_step=" << projStep << ",off_manifold_wait=" << offManifoldWait << ",collection_interval=" << collectInterval << ",nms=" << nMicroSteps << endl;
  fitted_coeffs_out << "n=" << n << ",proj_step=" << projStep << ",off_manifold_wait=" << offManifoldWait << ",collection_interval=" << collectInterval << ",nms=" << nMicroSteps << endl;
  pre_proj_degs_out.close();
  post_proj_degs_out.close();
  fitted_coeffs_out.close();

  
  if(init_type == "erdos") {
    initGraph();
  }
  else if(init_type == "complete") {
    init_complete_graph();
  }
  // initGraph() is default
  else {
    initGraph();
  }
    //set up save file and write header
    vector<double> forFile;
    forFile.push_back(projStep);
    forFile.push_back(offManifoldWait);
    forFile.push_back(nMicroSteps);
    forFile.push_back(collectInterval);
    vector<string> forFileStrs;
    forFileStrs.push_back("proj_step");
    forFileStrs.push_back("off_manifold_wait");
    forFileStrs.push_back("nms");
    forFileStrs.push_back("collection_interval");
    //open all the files
    // stringstream ss;
    // int folder_counter = 0;
    // string dir_base = "./csv_data";
    // string dir;
    // bool isdir = true;
    // struct stat stat_dir;
    // do {
    //   ss.str("");
    //   ss << dir_base << folder_counter << "/";
    //   folder_counter++;
    //   dir = ss.str();
    //   int check = stat(dir.c_str(), &stat_dir);
    //   if(check == -1) {
    // 	mkdir(dir.c_str(), 0700);
    // 	isdir = false;
    //   }
    // } while (isdir);
    // ofstream* paDataCPI = createFile("paDataCPI" + run_id, dir, forFile, forFileStrs);
    // ofstream* projData = createFile("projData" + run_id, dir, forFile, forFileStrs);
    // ofstream* eigVectData = createFile("eigVectData" + run_id, dir, forFile, forFileStrs);
    // ofstream* eigval_data = createFile("eigval_data" + run_id, dir, forFile, forFileStrs);
    // ofstream* deg_data = createFile("deg_data" + run_id, dir, forFile, forFileStrs);
    // ofstream* time_data = createFile("time_data" + run_id, dir , forFile, forFileStrs);
    //after waiting for the system to reach the slow manifold, collect data every collectInterval number of steps
    int totalSteps = 0;
    ofstream selfloops_out("./deg_cpi_data/selfloop_densities" + run_id + ".csv");
    selfloops_out << "n=" << n << ",proj_step=" << projStep << ",off_manifold_wait=" << offManifoldWait << ",collection_interval=" << collectInterval << ",nms=" << nMicroSteps << endl;
    while(totalSteps < nSteps) {
	int microStep;
	vector<graphData> toPlot;
	// vector<graphData> toProject;
	vector< vector<int> > degs_to_project;
	vector<double> time;
	vector<double> selfloop_densities;
	vector< vector<int> > degs_to_save;
	vector< int > times_to_save;

	// take one step and save this data, start microStep at 1
	// graphData* d = 
	step(true);
	// toPlot.push_back(*d);
	// degs_to_save.push_back(vector<int>(d->degSeq, d->degSeq+n));
	times_to_save.push_back(totalSteps);
	totalSteps++;
	for(microStep = 1; microStep < nMicroSteps; microStep++) {
	    if(microStep < offManifoldWait) {
		if((microStep+1)%save_interval == 0) {
		  // graphData* d = 
		  step(true);
		  // toPlot.push_back(*d);
		  // degs_to_save.push_back(vector<int>(d->degSeq, d->degSeq+n));
		  times_to_save.push_back(totalSteps);
		}
		else {
		    step(false);
		}
	    }
	    else {
	      int nOnManifoldSteps = microStep - offManifoldWait;
	      if((nOnManifoldSteps)%collectInterval == 0 || nOnManifoldSteps == nMicroSteps - offManifoldWait - 1) {
		if((microStep+1)%save_interval == 0) {
		  // graphData* d = 
		  step(true);
		  // toPlot.push_back(*d);
		  // degs_to_save.push_back(vector<int>(d->degSeq, d->degSeq+n));
		  times_to_save.push_back(totalSteps);
		  // toProject.push_back(*d);
		  time.push_back(totalSteps);
		  selfloop_densities.push_back(compute_selfloop_density());
		  degs_to_project.push_back(calcGraphProps::get_degrees(A, n));
		}
		else {
		  // toProject.push_back(*step(true));
		  time.push_back(totalSteps);
		  selfloop_densities.push_back(compute_selfloop_density());
		  degs_to_project.push_back(calcGraphProps::get_degrees(A, n));
		}
	      }
	      else {
		if((microStep+1)%save_interval == 0) {
		  // graphData* d =
		  step(true);
		  // toPlot.push_back(*d);
		  // degs_to_save.push_back(vector<int>(d->degSeq, d->degSeq+n));
		  times_to_save.push_back(totalSteps);
		}
		else {
		  step(false);
		}
	      }
	    }
	    totalSteps++;
	}

	// start MPI
	int size, rank;
	const int root = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	vector<int> new_degs(n);
	if(rank != root) {
	  int tag = 0;
	  for(vector< vector<int> >::iterator v = degs_to_project.begin(); v != degs_to_project.end(); v++) {
	    MPI_Send(&v->front(), n, MPI_INT, root, tag, MPI_COMM_WORLD);
	    tag++;
	  }
	}
	/* 
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	       probablty should use MPI_Gather
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	*/
	// gather all (use MPI_Gather()?) degs_to_project to
	// root process and average, then send out again to project with
	else {
	  for(int i = 0; i < degs_to_project.size(); i++) {
	    vector< vector<int> > other_degs(size-1, vector<int>(n));
	    int source = 1;
	    for(vector< vector<int> >::iterator v = other_degs.begin(); v != other_degs.end(); v++) {
	      MPI_Recv(&v->front(), n, MPI_INT, source, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	      source++;
	    }
	    // average across all processes
	    for(int j = 0; j < other_degs.size(); j++) {
	      for(int k = 0; k < n; k++) {
		degs_to_project[i][k] += other_degs[j][k];
	      }
	    }
	    for(int j = 0; j < n; j++) {
	      degs_to_project[i][j] /= size;
	    }
	  }
	  new_degs = project_degs(degs_to_project, time, projStep, run_id);
	}
	MPI_Bcast(&new_degs.front(), n, MPI_INT, root, MPI_COMM_WORLD);
	// end MPI
	  
	// vector<int> new_degs = project_degs(degs_to_project, time, projStep, run_id);
	init_graph_loosehh(new_degs);
	totalSteps += projStep;


	// save selfloop data
	for(int i = 0; i < selfloop_densities.size(); i++) {
	  selfloops_out << selfloop_densities[i] << endl;
	}


	// old projection method
	// //project collected data forward, save toPlot data, clear all vectors, update totalSteps
	// vector<vector<vector<double> > > data;
	// vector<vector<double> > mat;
	// vector<double> v;
	// for(vector<graphData>::iterator d = toProject.begin(); d != toProject.end(); d++) {
	//     data.push_back(mat);
	//     int i, j;
	//     for(i = 0; i < n; i++) {
	// 	data.back().push_back(v);
	// 	for(j = 0; j < n; j++) {
	// 	    (data.back()[i]).push_back((*d).A[i][j]);
	// 	}
	//     }
	// }
	// // project(data, time, projData, eigVectData, eigval_data);
	// int nPlotPts = toPlot.size();
	// //hooray for vlas
	// graphData toPlotAry[nPlotPts];
	// int i;
	// for(i = 0; i < nPlotPts; i++) {
	//     toPlotAry[i] = toPlot[i];
	// }
	// saveData(toPlotAry, nPlotPts, *paDataCPI);
	// save_degrees(degs_to_save, *deg_data);
	// saveData<int>(times_to_save, *time_data);
	// degs_to_save.clear();
	// times_to_save.clear();
	// toPlot.clear();

    }
    // *eigVectData << endl;
    // *eigval_data << endl;
    // paDataCPI->close();
    // projData->close();
    // eigVectData->close();
    // eigval_data->close();
    // delete paDataCPI;
    // delete projData;
    // delete eigVectData;
    // delete eigval_data;
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
void pamCPI::project(vector<vector<vector<double> > > &data, vector<double> &time, ofstream* projData, ofstream* eigVectData, ofstream* eigval_data) {
  //save data. this can be deleted after the damn thing runs properly
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
      line.push_back(50.0*i/n);
    }
    vector< vector< double > > coeffs_to_save;
    vector< vector< double > > eigval_tosave;
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
	eigval_tosave.push_back(leadingEigVect);
	for(j = 0; j < n; j++) {
	  delete[] eigVects[j];
	}
	delete[] eigVects;
	delete[] eigVals;


	// TESTING
	double recon_err = 0;
	double max_recon_err = 0;
	for(j = 0; j < n; j++) {
	  for(k = 0; k < n; k++) {
	    double err = abs(maxEigVal*leadingEigVect[j]*leadingEigVect[k] - currentAdjMat[j][k]);
	    recon_err += err;
	    if(err > max_recon_err) {
	      max_recon_err = err;
	    }
	  }
	}
	// cout << recon_err << " , " << max_recon_err << endl;


	sort(leadingEigVect.begin(), leadingEigVect.end());
	eigVectFittedCoeffs.push_back(fitCurves::fitFx(line, leadingEigVect, toFitEigVects));
	eigVectFittedCoeffs.back().push_back(maxEigVal);
	coeffs_to_save.push_back(eigVectFittedCoeffs.back());
	coeffs_to_save.back().push_back(time[i]);
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

    // TESTING
    // attempt to project only first coefficient (the offset) and eigenvalue
    // while averaging the others
    vector<double> to_average(nPts);
    for(i = 1; i < nCoeffs-1; i++) {
      for(j = 0; j < nPts; j++) {
	to_average[j] = (eigVectFittedCoeffs[j][i]);
      }
      newCoeffs[i] = utils::average(to_average);
    }
    // not pretty, but needed to make every time mod 100
    coeffs_to_save.back().back() = coeffs_to_save.back().back() + 1;

    coeffs_to_save.push_back(newCoeffs);
    coeffs_to_save.back().push_back(time.back() + 1 + projStep);
    double newEigVal = newCoeffs.back();
    newCoeffs.pop_back();
    nCoeffs = newCoeffs.size();
    vector<double> newEigVect;
    for(i = 0; i < n; i++) {
	double eval = 0;
	for(j = 0; j < nCoeffs; j++) {
	    eval += (*toFitEigVects[j])(line[i])*newCoeffs[j];
	}
	newEigVect.push_back(eval);
    }
    eigval_tosave.push_back(newEigVect);
    int **newA = new int*[n];
    for(i = 0; i < n; i++) {
	newA[i] = new int[n];
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
	}
    }

    initGraph(newA);

    vector<vector<double> > toSaveRecon;
    for(i = 0; i < n; i++) {
      toSaveRecon.push_back(vector<double>(n));
      for(j = 0; j < n; j++) {
	toSaveRecon[i][j] = A[i][j];
      }
    }
    for(i = 0; i < n; i++) {
	delete[] newA[i];
    }
    delete newA;
    //save data for visual comparison of reconstructions:
    vector<vector<double> > toSavePreRecon = data.back();
    saveData(toSavePreRecon, *projData);
    saveData(toSaveRecon, *projData);
    //save new eigvect
    sort(newEigVect.begin(), newEigVect.end());
    save_coeffs(coeffs_to_save, *eigVectData);
    save_coeffs(eigval_tosave, *eigval_data);
}

vector<int> pamCPI::project_degs(const std::vector< std::vector<int> >& deg_data, const std::vector<double>& times, const int proj_step, const string run_id) {
  // if proj_step is too small due negative degree restarts, simply return the current degree distribution
  // the choice of cutoff is arbitrary
  // if(proj_step < m*m/2) {
  //   return deg_data.back();
  // }

  const int n = deg_data[0].size();
  const int npts = deg_data.size();

  // add desired functions to fit to deg vs. index data
  fxs deg_fit_fns;
  fx const_offset = [] (double x) -> double { 
	return 1;
    };
    deg_fit_fns.push_back(const_offset);
    deg_fit_fns.push_back([] (double x) { return x;});
    deg_fit_fns.push_back([] (double x) { return x*x;});
    deg_fit_fns.push_back([] (double x) { return x*x*x;});
    deg_fit_fns.push_back([] (double x) { return x*x*x*x;});
    deg_fit_fns.push_back([] (double x) { return x*x*x*x*x;});
    
    vector<double> indices(n);
    for(int i = 0; i < n; i++) {
      indices[i] = i;
    }
    
  vector< vector<double> > deg_fitted_coeffs(npts, vector<double>());
    for(int i = 0; i < npts; i++) {
      vector<double> degrees(deg_data[i].begin(), deg_data[i].end());
      sort(degrees.begin(), degrees.end());
      deg_fitted_coeffs[i] = fitCurves::fitFx(indices, degrees, deg_fit_fns);
    }
    
    // fit each coefficient's evolution with a line
    const int ncoeffs = deg_fit_fns.size();
    fxs coeff_fit_fns;
    coeff_fit_fns.push_back(const_offset);
    coeff_fit_fns.push_back([] (double x) { return x;});
    vector< vector<double> > coeff_fitted_coeffs(ncoeffs);
    for(int i = 0; i < ncoeffs; i++) {
      vector<double> coeff_timecourse(npts);
      for(int j = 0; j < npts; j++) {
	coeff_timecourse[j] = deg_fitted_coeffs[j][i];
      }
      coeff_fitted_coeffs[i] = fitCurves::fitFx(times, coeff_timecourse, coeff_fit_fns);
    }
      
    // project data
    vector<double> new_coeffs(ncoeffs);
    double proj_time = times.back() + proj_step;
    for(int i = 0; i < ncoeffs; i++) {
      // evaluate each coefficient at the new time
      double eval = 0;
      for(int j = 0; j < coeff_fit_fns.size(); j++) {
	eval += (*coeff_fit_fns[j])(proj_time)*coeff_fitted_coeffs[i][j];
      }
      new_coeffs[i] = eval;
    }

    // reconstruct degree distribution
    vector<int> projected_degs(n);
    for(int i = 0; i < n; i++) {
      double eval = 0;
      for(int j = 0; j < ncoeffs; j++) {
	eval += (*deg_fit_fns[j])(indices[i])*new_coeffs[j];
      }
      projected_degs[i] = eval;
    }
    
    // rescale to approach approximately 'm' edges
    // allow maximum 0.1% disparity

    int degcount = 0;
    for(int i = 0; i < n; i++) {
      degcount += projected_degs[i];
    }
    double deg_scaling = 2.0*m/degcount;
    // while(std::abs(1 - deg_scaling) > 0.001) {
    for(int i = 0; i < n; i++) {
      projected_degs[i] = deg_scaling*projected_degs[i];
    }
    //   degcount = 0;
    //   for(int i = 0; i < n; i++) {
    // 	degcount += projected_degs[i];
    //   }
    //   deg_scaling = 2.0*m/degcount;
    // }

    // ensure all degrees are non-negative
    // if a negative value is found, re-try with half the projection step
    for(int i = 0; i < n; i++) {
      if(projected_degs[i] < 0) {

	// TESTING
	projected_degs[i] = 0;

	// return project_degs(deg_data, time, proj_step/2);
      }
    }

   // TESTING
    degcount = 0;
    for(int i = 0; i < n; i++) {
      degcount += projected_degs[i];
    }

    // cout << "degree difference: " << degcount - 2*m << endl;

    ofstream times_out("./deg_cpi_data/times" + run_id + ".csv", ios_base::app);
    ofstream pre_proj_degs_out("./deg_cpi_data/pre_proj_degs" + run_id + ".csv", ios_base::app);
    ofstream post_proj_degs_out("./deg_cpi_data/post_proj_degs" + run_id + ".csv", ios_base::app);
    ofstream fitted_coeffs_out("./deg_cpi_data/fitted_coeffs" + run_id + ".csv", ios_base::app);
    // truly abhorrent initialization
    saveData(times, times_out);
    save_coeffs(vector< vector<double> >(1, std::vector<double>(deg_data.back().begin(), deg_data.back().end())), pre_proj_degs_out);
    save_coeffs(vector< vector<double> >(1, std::vector<double>(projected_degs.begin(), projected_degs.end())), post_proj_degs_out);
    // save_coeffs(vector< vector<double> >(1, deg_fitted_coeffs.back()), fitted_coeffs_out);
    save_coeffs(deg_fitted_coeffs, fitted_coeffs_out);
    return projected_degs;
}
