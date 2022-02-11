/*
 * Copyright 319 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */

 #include <iostream>

 #include <sstream>
 #include <fstream>
 #include <vector>
 #include <string>
 #include <iomanip>
 #include <cstring>
 #include <cmath>
 #include <cstdlib>
// standard
#include <stdio.h>
#include <stdlib.h>
// acados
#include "acados/utils/print.h"
#include "acados/utils/math.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"
#include "acados_solver_mhe_pendulum_ode.h"

using namespace std;

bool readDataFromFile( const char* fileName, vector< vector< double > >& data )
{
	ifstream file( fileName );
	string line;

	if ( file.is_open() )
	{
		while( getline(file, line) )
		{
			istringstream linestream( line );
			vector< double > linedata;
			double number;

			while( linestream >> number )
			{
				linedata.push_back( number );
			}

			data.push_back( linedata );
		}

		file.close();
	}
	else
		return false;

	return true;
}

int main()
{
  vector< vector< double > > simX;
  if (readDataFromFile("../solver/mhe/simX.txt", simX) == false)
  {
    cout << "Cannot read simX" << endl;
    return EXIT_FAILURE;
  }

  vector< vector< double > > simU;
  if (readDataFromFile("../solver/mhe/simU.txt", simU) == false)
  {
    cout << "Cannot read simU" << endl;
    return EXIT_FAILURE;
  }

  vector< vector< double > > simY;
  if (readDataFromFile("../solver/mhe/simY.txt", simY) == false)
  {
    cout << "Cannot read simY" << endl;
    return EXIT_FAILURE;
  }

  printf("\n--- simX ---\n");
  for(int i=0; i<simX.size(); ++i) {
    for(int j=0; j<simX[i].size(); ++j) {
      std::cout << simX[i][j] << ' ';
    }
    std::cout << std::endl;
  }

  printf("\n--- simU ---\n");
  for(int i=0; i<simU.size(); ++i)
  std::cout << simU[i][0] << ' ';

  printf("\n--- simY ---\n");
  for(int i=0; i<simY.size(); ++i) {
    for(int j=0; j<simY[i].size(); ++j) {
      std::cout << simY[i][j] << ' ';
    }
    std::cout << std::endl;
  }

    int status = 0;
    status = acados_create();

    if (status)
    {
        printf("acados_create() returned status %d. Exiting.\n", status);
        exit(1);
    }

    // initial condition
    int idxbx0[0];

    double lbx0[0];
    double ubx0[0];

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxbx", idxbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx", lbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx", ubx0);

    double yref[12] = {0};
    // initialization for state values
    double x_init[4];
    x_init[0] = 0.0;
    x_init[1] = 0.0;
    x_init[2] = 0.0;
    x_init[3] = 0.0;

    // initial value for control input
    double u0[4];
    u0[0] = 0.0;
    u0[1] = 0.0;
    u0[2] = 0.0;
    u0[3] = 0.0;

    // arrival cost mean
    double x0_bar[4] = {0.0, 3.141592653589793, 0.0, 0.0};

    // set parameters
    double p[1];

    p[0] = 0;

    for (int ii = 0; ii < 3; ii++)
    {
        expl_ode_fun[ii].set_param(expl_ode_fun+ii, p);
        forw_vde_casadi[ii].set_param(forw_vde_casadi+ii, p);
    }
    for (int ii = 0; ii < 3; ii++) {
    }


    // prepare evaluation
    int NTIMINGS = 1;
    double min_time = 1e12;
    double kkt_norm_inf;
    double elapsed_time;
    int sqp_iter;

    double xtraj[4 * (3+1)];
    double utraj[4 * (3)];


    // solve ocp in loop
    int rti_phase = 0;

    // // set arrival cost penalty (W at stage 0)
    // vector< vector< double > > W0;
    // if (readDataFromFile("../solver/mhe/W0.txt", W0) == false)
    // {
    //     cout << "Cannot read W0" << endl;
    //     return EXIT_FAILURE;
    // }
        /**** Cost ****/
    #define NY 12
    double W[NY*NY];
    
    W[0+(NY) * 0] = 0.1;
    W[0+(NY) * 1] = 0;
    W[0+(NY) * 2] = 0;
    W[0+(NY) * 3] = 0;
    W[0+(NY) * 4] = 0;
    W[0+(NY) * 5] = 0;
    W[0+(NY) * 6] = 0;
    W[0+(NY) * 7] = 0;
    W[0+(NY) * 8] = 0;
    W[0+(NY) * 9] = 0;
    W[0+(NY) * 10] = 0;
    W[0+(NY) * 11] = 0;
    W[1+(NY) * 0] = 0;
    W[1+(NY) * 1] = 0.1;
    W[1+(NY) * 2] = 0;
    W[1+(NY) * 3] = 0;
    W[1+(NY) * 4] = 0;
    W[1+(NY) * 5] = 0;
    W[1+(NY) * 6] = 0;
    W[1+(NY) * 7] = 0;
    W[1+(NY) * 8] = 0;
    W[1+(NY) * 9] = 0;
    W[1+(NY) * 10] = 0;
    W[1+(NY) * 11] = 0;
    W[2+(NY) * 0] = 0;
    W[2+(NY) * 1] = 0;
    W[2+(NY) * 2] = 0.1;
    W[2+(NY) * 3] = 0;
    W[2+(NY) * 4] = 0;
    W[2+(NY) * 5] = 0;
    W[2+(NY) * 6] = 0;
    W[2+(NY) * 7] = 0;
    W[2+(NY) * 8] = 0;
    W[2+(NY) * 9] = 0;
    W[2+(NY) * 10] = 0;
    W[2+(NY) * 11] = 0;
    W[3+(NY) * 0] = 0;
    W[3+(NY) * 1] = 0;
    W[3+(NY) * 2] = 0;
    W[3+(NY) * 3] = 0.1;
    W[3+(NY) * 4] = 0;
    W[3+(NY) * 5] = 0;
    W[3+(NY) * 6] = 0;
    W[3+(NY) * 7] = 0;
    W[3+(NY) * 8] = 0;
    W[3+(NY) * 9] = 0;
    W[3+(NY) * 10] = 0;
    W[3+(NY) * 11] = 0;
    W[4+(NY) * 0] = 0;
    W[4+(NY) * 1] = 0;
    W[4+(NY) * 2] = 0;
    W[4+(NY) * 3] = 0;
    W[4+(NY) * 4] = 0.1;
    W[4+(NY) * 5] = 0;
    W[4+(NY) * 6] = 0;
    W[4+(NY) * 7] = 0;
    W[4+(NY) * 8] = 0;
    W[4+(NY) * 9] = 0;
    W[4+(NY) * 10] = 0;
    W[4+(NY) * 11] = 0;
    W[5+(NY) * 0] = 0;
    W[5+(NY) * 1] = 0;
    W[5+(NY) * 2] = 0;
    W[5+(NY) * 3] = 0;
    W[5+(NY) * 4] = 0;
    W[5+(NY) * 5] = 0.1;
    W[5+(NY) * 6] = 0;
    W[5+(NY) * 7] = 0;
    W[5+(NY) * 8] = 0;
    W[5+(NY) * 9] = 0;
    W[5+(NY) * 10] = 0;
    W[5+(NY) * 11] = 0;
    W[6+(NY) * 0] = 0;
    W[6+(NY) * 1] = 0;
    W[6+(NY) * 2] = 0;
    W[6+(NY) * 3] = 0;
    W[6+(NY) * 4] = 0;
    W[6+(NY) * 5] = 0;
    W[6+(NY) * 6] = 0.1;
    W[6+(NY) * 7] = 0;
    W[6+(NY) * 8] = 0;
    W[6+(NY) * 9] = 0;
    W[6+(NY) * 10] = 0;
    W[6+(NY) * 11] = 0;
    W[7+(NY) * 0] = 0;
    W[7+(NY) * 1] = 0;
    W[7+(NY) * 2] = 0;
    W[7+(NY) * 3] = 0;
    W[7+(NY) * 4] = 0;
    W[7+(NY) * 5] = 0;
    W[7+(NY) * 6] = 0;
    W[7+(NY) * 7] = 0.1;
    W[7+(NY) * 8] = 0;
    W[7+(NY) * 9] = 0;
    W[7+(NY) * 10] = 0;
    W[7+(NY) * 11] = 0;
    W[8+(NY) * 0] = 0;
    W[8+(NY) * 1] = 0;
    W[8+(NY) * 2] = 0;
    W[8+(NY) * 3] = 0;
    W[8+(NY) * 4] = 0;
    W[8+(NY) * 5] = 0;
    W[8+(NY) * 6] = 0;
    W[8+(NY) * 7] = 0;
    W[8+(NY) * 8] = 100;
    W[8+(NY) * 9] = 0;
    W[8+(NY) * 10] = 0;
    W[8+(NY) * 11] = 0;
    W[9+(NY) * 0] = 0;
    W[9+(NY) * 1] = 0;
    W[9+(NY) * 2] = 0;
    W[9+(NY) * 3] = 0;
    W[9+(NY) * 4] = 0;
    W[9+(NY) * 5] = 0;
    W[9+(NY) * 6] = 0;
    W[9+(NY) * 7] = 0;
    W[9+(NY) * 8] = 0;
    W[9+(NY) * 9] = 100;
    W[9+(NY) * 10] = 0;
    W[9+(NY) * 11] = 0;
    W[10+(NY) * 0] = 0;
    W[10+(NY) * 1] = 0;
    W[10+(NY) * 2] = 0;
    W[10+(NY) * 3] = 0;
    W[10+(NY) * 4] = 0;
    W[10+(NY) * 5] = 0;
    W[10+(NY) * 6] = 0;
    W[10+(NY) * 7] = 0;
    W[10+(NY) * 8] = 0;
    W[10+(NY) * 9] = 0;
    W[10+(NY) * 10] = 100;
    W[10+(NY) * 11] = 0;
    W[11+(NY) * 0] = 0;
    W[11+(NY) * 1] = 0;
    W[11+(NY) * 2] = 0;
    W[11+(NY) * 3] = 0;
    W[11+(NY) * 4] = 0;
    W[11+(NY) * 5] = 0;
    W[11+(NY) * 6] = 0;
    W[11+(NY) * 7] = 0;
    W[11+(NY) * 8] = 0;
    W[11+(NY) * 9] = 0;
    W[11+(NY) * 10] = 0;
    W[11+(NY) * 11] = 100;
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, 0, "W", W);


    for (int ii = 0; ii < NTIMINGS; ii++)
    {
        printf("nlp N = %d\n", nlp_dims->N);
        // initialize solution
        // TODO: loop to N with < or <=?!
        for (int i = 0; i < nlp_dims->N; i++)
        {
            yref[0] = simY[i][0];
            yref[1] = simY[i][1];
            yref[2] = simY[i][2];
            yref[3] = simY[i][3];

            yref[8] = x0_bar[0];
            yref[9] = x0_bar[1];
            yref[10] = x0_bar[2];
            yref[11] = x0_bar[3];

            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "x", x_init);
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "u", u0);
            ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "yref", yref);
            printf("uval, we wanna set: %e, at stage %d\n", simU[i][0], i);
            acados_update_params(i, &simU[i][0], 1);
            printf("param update succes\n");
        }
        
        ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "rti_phase", &rti_phase);
        status = acados_solve();
        printf("\nIN main, after solve\n");
        ocp_nlp_get(nlp_config, nlp_solver, "time_tot", &elapsed_time);
        min_time = MIN(elapsed_time, min_time);
    }

    /* print solution and statistics */
    for (int ii = 0; ii <= nlp_dims->N; ii++)
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, ii, "x", &xtraj[ii*4]);
    for (int ii = 0; ii < nlp_dims->N; ii++)
        ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, ii, "u", &utraj[ii*4]);

    printf("\n--- xtraj ---\n");
    d_print_exp_tran_mat( 4, 3+1, xtraj, 4 );
    printf("\n--- utraj ---\n");
    d_print_exp_tran_mat( 4, 3, utraj, 4 );
    // ocp_nlp_out_print(nlp_solver->dims, nlp_out);

    printf("\nsolved ocp %d times, solution printed above\n\n", NTIMINGS);

    if (status == ACADOS_SUCCESS)
    {
        printf("acados_solve(): SUCCESS!\n");
    }
    else
    {
        printf("acados_solve() failed with status %d.\n", status);
    }

    // get solution
    ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 0, "kkt_norm_inf", &kkt_norm_inf);
    ocp_nlp_get(nlp_config, nlp_solver, "sqp_iter", &sqp_iter);

    acados_print_stats();

    printf("\nSolver info:\n");
    printf(" SQP iterations %2d\n minimum time for %d solve %f [ms]\n KKT %e\n",
           sqp_iter, NTIMINGS, min_time*1000, kkt_norm_inf);

    
    // save xtarj to a file
    FILE *f = fopen("simXest.txt", "w");
    for (int i=0; i<3+1; ++i) {
        for (int j=0; j<4; ++j) {
            fprintf(f, "%.18e ", xtraj[i * 4 + j]);
        }
        fprintf(f, "\n");
    }
    for (int i=0; i<sizeof(xtraj); ++i)
    {
        
    }
//     fwrite(xtraj, sizeof(xtraj), 1, f);
    fclose(f);

    // free solver
    status = acados_free();
    if (status) {
        printf("acados_free() returned status %d. \n", status);
    }

    return status;
}
