#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <float.h>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>

#include <fstream>
#include <sstream>
#include <iostream>

#ifdef DEBUG
   #ifdef NDEBUG
      #undef NDEBUG
   #endif
#else
   #define NDEBUG
#endif
#include <assert.h>

#include <clock.h>

#include <cklib.h>
#include <rk.h>
#include <ros.h>
#include <sdirk.h>

#define runtime_assert( cmp ) { \
   if (! (cmp) ) { \
      fprintf(stderr,"Run-time assertion: %s:%d\n\t%s\n", __FILE__, __LINE__, #cmp); return 1; } \
   }


#ifdef _OPENMP
# include <omp.h>
#endif

#include <Vector.h>
#include <ck.hpp>

#define getTempIndex(neq) ( (neq)-1 )
#define getFirstSpeciesIndex(neq) ( 0 )

#define ENABLE_SIMD
#ifdef ENABLE_SIMD
# include "simd_dtypes.hpp"
# include "simd_cklib.hpp"
# include "simd_solvers.hpp"
#endif

void ckrhs (const double p,
            const double T,
            const double y[],
                  double *Tdot,
                  double ydot[],
            const ckdata_t *ck,
                  double rwk[])
{
   const int kk = ck->n_species;

   /* Compute molar reaction rate. (mol / (s*cm^3) */
   ckwyp (p, T, y, ydot, ck, rwk);

   /* Compute mixture Cp (ergs / gm*K) */
   const double cp_mix = ckcpbs(T, y, ck);

   /* Compute species enthalpy (ergs / K) */
   ckhms(T, rwk, ck);

   /* Extract the molecular weights of the species ... this could just be a pointer. */
   const double rho = ckrhoy(p, T, y, ck);

   *Tdot = 0.0;
   for (int k = 0; k < kk; ++k)
   {
      /* Convert from molar to mass units. */
      ydot[k] *= ck->sp_mwt[k];
      ydot[k] /= rho;

      /* Sum up the net enthalpy change. */
      *Tdot -= (rwk[k] * ydot[k]);
   }

   *Tdot /= cp_mix;

   return;
}

struct cklib_functor
{
   const ckdata_t *m_ckptr;

   double m_pres; // dynes/cm^2

   int m_lenrwk;
   VectorType<double,64> m_rwk;

   cklib_functor(const ckdata_t *ckptr, const double constant_p = 1.013250e6 ) :
      m_ckptr(ckptr),
      m_pres(constant_p),
      m_lenrwk(0), m_rwk()
   {
      m_lenrwk = ck_lenrwk(m_ckptr);
      m_rwk.resize( m_lenrwk );
   }

   ~cklib_functor()
   {
   }

   double getPressure(void) const { return this->m_pres; }

   void operator() (const int &neq, double y[], double f[])
   {
      this->operator()(neq, 0.0, y, f);
   }
   void operator() (const int &neq, const double &time, double y[], double f[])
   {
      this->rhs(neq, time, y, f);
   }
   int rhs (const int &neq, const double &time, double y[], double f[])
   {
      const int kk = this->m_ckptr->n_species;

      const double T = y[ getTempIndex(neq) ];
      const double *yk = y + getFirstSpeciesIndex(neq);

      double *ykdot = f + getFirstSpeciesIndex(neq);
      double *Tdot = f + getTempIndex(neq);

      ckrhs ( this->m_pres, T, yk, Tdot, ykdot, this->m_ckptr, this->m_rwk.getPointer() );

      return 0;
   }

   int jac (const int &neq, double t, double y[], double fy[])
   {
      runtime_assert(true);
   }
};

static int rhs_func (const int neq, const double time, double y[], double f[], void *v_user)
{
   cklib_functor *ck = (cklib_functor *) v_user;
   return ck->rhs( neq, time, y, f );
}
int cklib_callback (const int neq, const double time, double y[], double f[], void *v_user)
{
   fprintf(stderr,"Still calling cklib_callback\n");
   exit(-1);
   return -1;
}
static int jac_func (const int neq, const double time, double y[], double fy[], void *v_user)
{
   cklib_functor *ck = (cklib_functor *) v_user;
   return ck->jac( neq, time, y, fy );
}

#ifdef USE_TCHEM
#warning 'Enabling TChem'
namespace TChem
{

extern "C"
{
# include <TC_interface.h>
# include <TC_defs.h>
}

struct Functor
{
   Functor( const double constant_p = 1.013250e6 )
   {
      std::string cheminp = "chem.inp";
      std::string thermdat = "therm.dat";
      TC_initChem (const_cast<char*>(cheminp.c_str()), const_cast<char*>(thermdat.c_str()), 0, 0.2);

      /* Internal TChem variables */
      TC_abstol = 1.0e-9;
      TC_reltol = 1.0e-11;

      TC_setThermoPres (constant_p/10.);

      printf("TChem initialized with pressure= %f %d %d\n", constant_p, TC_getNspec(), TC_getNreac());
   }

   void rhs (const int neq, const double t, double *y, double *fy)
   {
#if ( getTempIndex(0) == 0 )
# warning 'Using TC_getSrc directly'
      TC_getSrc ( y, neq, fy );
#else
      std::vector<double> ytmp(neq), fytmp(neq);

      const int nsp = neq-1;

      ytmp[0] = y[nsp];
      for (int i = 0; i < nsp; ++i)
      {
         ytmp[i+1] = y[i];
      }

      TC_getSrc ( &ytmp[0], neq, &fytmp[0] );

      fy[nsp] = fytmp[0];
      for (int i = 0; i < nsp; ++i)
         fy[i] = fytmp[i+1];
#endif
   }
   int rhs (const int neq, const double t, double *y, double *fy, void *vdata)
   {
      this->rhs (neq, t, y, fy);
      return 0;
   }

   int operator() (const int neq, const double t, double *y, double *fy)
   {
      this->rhs (neq, t, y, fy);
      return 0;
   }

   int jac (const int &neq, double t, double y[], double fy[])
   {
      runtime_assert(true);
   }

};

} // namespace TChem
#endif

#ifdef USE_PYJAC
#warning 'Enabling pyJac'
namespace PyJac
{

extern "C"
{
# include <chem_utils.h>
# include <dydt.h>
# include <header.h>
# include <jacob.h>
# include <mass_mole.h>
# include <mechanism.h>
# include <rates.h>
# include <sparse_multiplier.h>
}

struct Functor
{
   double pres;

   Functor( const double constant_p = 1.013250e6 )
      : pres(constant_p/10.0)
   {
      printf("PyJac initialized with pressure= %f %d\n", constant_p, (NSP));
   }

   void rhs (const int neq, const double t, double *y, double *fy)
   {
#if ( getTempIndex(0) == 0 )
# warning 'Using PyJac directly (ignoring N v. (N-1) form.'
      dydt ( t, this->pres, y, fy );
#else
      std::vector<double> ytmp(neq), fytmp(neq);

      const int nsp = neq-1;

      ytmp[0] = y[nsp];
      for (int i = 0; i < nsp; ++i)
      {
         ytmp[i+1] = y[i];
      }

      dydt ( t, this->pres, &ytmp[0], &fytmp[0] );

      fy[nsp] = fytmp[0];
      for (int i = 0; i < nsp; ++i)
         fy[i] = fytmp[i+1];
#endif
   }
   int rhs (const int neq, const double t, double *y, double *fy, void *vdata)
   {
      this->rhs (neq, t, y, fy);
      return 0;
   }

   int operator() (const int neq, const double t, double *y, double *fy)
   {
      this->rhs (neq, t, y, fy);
      return 0;
   }

   int jac (const int &neq, double t, double y[], double fy[])
   {
      runtime_assert(true);
   }
};

} // namespace PyJac
#endif

#ifdef USE_SUNDIALS
# include <cv_integrator.h>

template <typename ValueType, typename ScalarType, class Solver, class Func>
std::tuple<int,int,int,int>
     cv_driver (const int         neq,
                      ValueType   u_in[],
                const ScalarType& t_stop,
                      Solver&     solver,
                      Func&       func,
                const bool write_data = false)
{
   double time_start = WallClock();

   ValueType *u = u_in;

   const int nsteps = (write_data == true) ? 1000 : 20;
   ScalarType t0 = 0.0;
   ScalarType dt = t_stop / double(nsteps);

   ValueType t(t0);

   int nst = 0;
   int nni = 0;
   int nlu = 0;
   int nfe = 0;
   int nje = 0;

   int ierr = solver.init (t, t_stop, u, func);
   if (ierr)
   {
      std::cerr << "Error in CV::Integrator::init()" << std::endl;
      return std::make_tuple( nst, nfe, nje, nni );
   }

   for (int i = 0; i < nsteps; ++i)
   {
      ScalarType t_next = (i == nsteps-1) ? t_stop : fmin(t + dt, t_stop);
      double _t0 = WallClock();

      const int itask = CV_NORMAL;
      solver.solve(t, t_next, u, func, itask);

      if (nsteps > 1)
      {
         double _t1 = WallClock();
         //std::cout << "step: " << i << std::endl;
         //std::cout << "   t: " << t << std::endl;
         //std::cout << "   T: " << u[neq-1] << std::endl;
         //std::cout << " nst: " << solver.nst << std::endl;
         //std::cout << "time: " << 1000*(_t1-_t0) << std::endl;
         printf("%5d %g, %g, %d, %d, %g\n", i, t, u[neq-1], solver.nst, nst, 1000*(_t1-_t0));
      }

      nst += solver.nst;
      nfe += solver.nfe;
      nje += solver.nje;
      nni += solver.nni;
      nlu += solver.nlu;
   }

   double time_stop = WallClock();

   return std::make_tuple( nst, nfe, nje, nni );
}

#endif

int main (int argc, char* argv[])
{
   std::string ckname("ck.bin");

   typedef enum { None, RK, ROS, SDIRK, CV } SolverTags;

   int numProblems = 1;
   SolverTags solverTag = None;
   double t_stop = 0.001;
   char *load_from_file = NULL;
   bool randomize_input = false;
   bool enableTChem = false;
   bool enablePyJac = false;

   for (int index = 1; index < argc; ++index)
   {
      printf("argv[%d] = %s\n", index, argv[index]);

      if (strcmp(argv[index], "-ck") == 0)
      {
         index++;
         runtime_assert( index < argc );
         ckname = std::string(argv[index]);
      }
      else if (strcmp(argv[index], "-np") == 0)
      {
         index++;
         runtime_assert (index < argc);
         runtime_assert( isdigit(*argv[index]) ) ;
         numProblems = atoi(argv[index]);
      }
      else if (strcmp(argv[index], "-tstop") == 0)
      {
         index++;
         runtime_assert( index < argc );
         runtime_assert( isdigit(*argv[index]) ) ;
         t_stop = atof(argv[index]);
      }
      else if (strcmp(argv[index], "-read-from") == 0)
      {
         index++;
         runtime_assert( index < argc );
         load_from_file = argv[index];
      }
      else if (strcmp(argv[index], "-tchem") == 0)
      {
         enableTChem = true;
      }
      else if (strcmp(argv[index], "-pyjac") == 0)
      {
         enablePyJac = true;
      }
      else if (strcmp(argv[index], "-ros") == 0)
      {
         solverTag = ROS;
      }
      else if (strcmp(argv[index], "-rk") == 0)
      {
         solverTag = RK;
      }
      else if (strcmp(argv[index], "-sdirk") == 0)
      {
         solverTag = SDIRK;
      }
      else if (strcmp(argv[index], "-cv") == 0)
      {
         solverTag = CV;
      }
      else if (strcmp(argv[index], "-rand") == 0)
      {
         randomize_input = true;
      }
   }

   FILE *ckfile = fopen(ckname.c_str(),"r");
   if (not(ckfile))
   {
      fprintf(stderr,"error opening ckfile %s", ckname.c_str());
      return 1;
   }

   CK::CKData ck;

   CK::ckinit (ck, ckfile);

   fclose(ckfile);

   const double p = CK::PA;

   ckdata_t *ckptr = NULL;

   // Convert C++ CK object to c99 version.
   {
      std::vector<char *> sp_name_(ck.n_species);
      for (int k = 0; k < ck.n_species; ++k)
         sp_name_[k] = const_cast<char*>(ck.sp_name[k].c_str());

      std::vector<int> rx_falloff_type_(ck.n_falloff+1);
      for (int i = 0; i < ck.n_falloff; ++i)
      {
         int k = ck.rx_falloff_idx[i];
         int type = 0;
         if (!(ck.rx_info[k] & CK::rx_flag_falloff))
            printf("not a falloff %d %d\n", i, k);
         if (ck.rx_info[k] & CK::rx_flag_falloff_sri)
         {
            type = 1;
            if (ck.rx_info[k] & CK::rx_flag_falloff_sri5)
               type = 2;
         }
         else if (ck.rx_info[k] & CK::rx_flag_falloff_troe)
         {
            type = 3;
            if (ck.rx_info[k] & CK::rx_flag_falloff_troe4)
               type = 4;
         }
         rx_falloff_type_[i] = type;
         //printf("%d %d %d %d\n", i, k, ck.rx_info[k], type);
      }

      ckptr = ck_create (
         ck.n_species,
         &sp_name_[0], ck.sp_mwt.getPointer(), ck.th_tmid.getPointer(), ck.th_alo.getPointer(), ck.th_ahi.getPointer(),
         ck.n_reactions,
         ck.rx_A.getPointer(), ck.rx_b.getPointer(), ck.rx_E.getPointer(), ck.rx_nu.getPointer(), ck.rx_nuk.getPointer(),
         ck.n_reversible_reactions,
         ck.rx_rev_idx.getPointer(), ck.rx_rev_A.getPointer(), ck.rx_rev_b.getPointer(), ck.rx_rev_E.getPointer(),
         ck.n_irreversible_reactions,
         ck.rx_irrev_idx.getPointer(),
         ck.n_thdbdy,
         ck.rx_thdbdy_idx.getPointer(), ck.rx_thdbdy_offset.getPointer(), ck.rx_thdbdy_spidx.getPointer(), ck.rx_thdbdy_alpha.getPointer(),
         ck.n_falloff,
         ck.rx_falloff_idx.getPointer(), &rx_falloff_type_[0], ck.rx_falloff_spidx.getPointer(), ck.rx_falloff_params.getPointer());
   }

   const int kk = ck.n_species;
   const int neq = kk+1;

   std::cout << "numProblems (requested) = " << numProblems << std::endl;

   std::vector<double> u_in;

   // Load problem data from a file.
   if ( load_from_file )
   {
      int np = -1;

      // First, load all of the file into a list of strings for each line.
      std::vector <std::vector <std::string> > file_data;
      std::ifstream infile( load_from_file );

      while (infile)
      {
         std::string s;
         if (!std::getline( infile, s )) break;

         //std::cout << "line: " << s << std::endl;

         std::istringstream ss( s );
         std::vector <std::string> record;

         while (ss)
         {
            std::string s;
            if (!std::getline( ss, s, ',' )) break;
            record.push_back( s );
            //std::cout << s<< std::endl;
            //size_t pos = 0;
            //while (s[pos] == ' ') pos++;
            //record.push_back( std::string(s, pos, std::string::npos) );
         }

         file_data.push_back( record );
      }
      if (!infile.eof())
      {
         std::cerr << "Fooey!\n";
         return -1;
      }

      np = atoi( file_data[0][0].c_str() );
      if ( np > file_data.size()-1 )
      {
         fprintf(stderr,"Input profile file error: fewer lines than specified %d %d\n", np, file_data.size()-1);
         return -1;
      }

      int nsp = atoi( file_data[0][1].c_str() );
      if (nsp != ck.n_species)
      {
         fprintf(stderr,"Input profile file error: # species not correct %d %d\n", nsp, ck.n_species);
         return -1;
      }

      printf("numProblems (in) = %d %d\n", numProblems, np);

      std::vector< int   > order(np);

      for (int i = 0; i < np; ++i)
         order[i] = i;

      if ( randomize_input )
      {
         std::cout << "Randomizing input order\n";

         std::vector< float > f(np);

         auto comp = [&](const int &a, const int &b)
            {
               return f[ a ] < f[ b ] ;
            };

         for (int i = 0; i < np; ++i)
            f[i] = float( rand() ) / RAND_MAX;

         std::sort( order.begin(), order.end(), comp );

         //for (int i = 0; i < np; ++i)
         //   printf("%d %d %f\n", i, order[i], f[order[i]]);
      }

      //if (np < numProblems)
      //   numProblems = np;

      //t_stop = dt;

      u_in.resize( neq*numProblems );

      for (int i = 0; i < std::min(np, numProblems); ++i)
      {
         const int line_index = order[i];

         double x_, v_ = 0, T_, p_;
         std::vector<double> yk_(ck.n_species);

         x_ = strtod( file_data[line_index+1][0].c_str(), NULL );
         T_ = strtod( file_data[line_index+1][1].c_str(), NULL );
         p_ = strtod( file_data[line_index+1][2].c_str(), NULL );

         for (int k = 0; k < ck.n_species; ++k)
            yk_[k] = strtod( file_data[line_index+1][k+3].c_str(), NULL );

         const int offset = i * neq;
         u_in[ offset + getTempIndex(neq) ] = T_;
         for (int k = 0; k < ck.n_species; ++k)
            u_in[ offset + getFirstSpeciesIndex(neq) + k ] = yk_[k];
      }

      /* Now, fill in the rest of the slots */
      if (np < numProblems)
      {
         printf("replicating inputs to fill %d problems\n", numProblems);
         for (int i = np; i < numProblems; ++i)
         {
            int j = i % np;
            memcpy( &u_in[i*neq], &u_in[j*neq], sizeof(double)*neq);
         }
      }
   }
   else
   {
      // Define a simple ignition problem with the a common fuel.
      int iH2=-1, iO2=-1, iN2=-1, iCH4=-1;
      for (int k = 0; k < ck.n_species; ++k)
         if      (ck.sp_name[k].compare("O2")  == 0) iO2  = k;
         else if (ck.sp_name[k].compare("H2")  == 0) iH2  = k;
         else if (ck.sp_name[k].compare("N2")  == 0) iN2  = k;
         else if (ck.sp_name[k].compare("CH4") == 0) iCH4 = k;

      printf("iH2=%d, iO2=%d, iN2=%d, iCH4=%d\n", iH2, iO2, iN2, iCH4);

      std::vector<double> x(ck.n_species);
      std::vector<double> y(ck.n_species);

      x[iH2] = 2.0;
      x[iO2] = 1.0;
      if (iN2 != -1) x[iN2] = 4.0;

      const double T = 1000.0; // K
      double x_sum(0);
      for (int k = 0; k < ck.n_species; ++k)
         x_sum += x[k];

      for (int k = 0; k < ck.n_species; ++k)
         x[k] /= x_sum;

      // Mole to mass fractions.
      CK::ckxty (&x[0], &y[0], ck);
      const double rho = CK::ckrhoy(p, T, &y[0], ck);
      printf("cp_mass=%e, cv_mass=%e, h_mass=%e, u_mass=%e, rho=%e, H2=%f, O2=%f, N2=%f, T=%f\n", CK::ckcpbs(T,&y[0],ck), CK::ckcvbs(T,&y[0],ck), CK::ckhbms(T,&y[0],ck), CK::ckubms(T,&y[0],ck), rho, x[iH2], x[iO2], (iN2 != -1) ? x[iN2] : -1.0, T);

      numProblems = 1;
      u_in.resize( neq );
      u_in[ getTempIndex(neq) ] = T;
      for (int k = 0; k < ck.n_species; ++k)
         u_in[ getFirstSpeciesIndex(neq) + k] = y[k];
   }

   VectorType<double> cv_err(neq), cv_ref(neq);
   for (int i = 0; i < neq; ++i)
      cv_err[i] = cv_ref[i] = 0;

   VectorType<double> u(neq);

   cklib_functor cklib_func( ckptr, p );

#ifdef ENABLE_SIMD
   {
      SIMD::test_simd_rhs( numProblems, &u_in[0], cklib_func, ckptr );
      if (solverTag == RK) {
         SIMD::simd_rk_driver( numProblems, &u_in[0], t_stop, cklib_func, rhs_func, ckptr );
         return 1;
      }
      else if (solverTag == ROS) {
         SIMD::simd_ros_driver( numProblems, &u_in[0], t_stop, cklib_func, rhs_func, ckptr );
         return 1;
      }
      else if (solverTag == SDIRK) {
         SIMD::simd_sdirk_driver( numProblems, &u_in[0], t_stop, cklib_func, rhs_func, ckptr );
         return 1;
      }
   }
#endif

#ifdef USE_SUNDIALS
   CV::Integrator< cklib_functor > cv_obj( neq );
#endif

#ifdef USE_TCHEM
   CV::Integrator< TChem::Functor > cv_tchem( neq );
   TChem::Functor tchem_func;
   if (enableTChem)
   {
      VectorType<double> tc_out(neq), ck_out(neq);

      cklib_func( neq, 0.0, &u_in[0], ck_out.getPointer());
      tchem_func( neq, 0.0, &u_in[0], tc_out.getPointer());
      for (int i = 0; i < neq; ++i)
      {
         double refval = std::fabs(tc_out[i]);
         if ( refval <= 1e-20 ) refval = 1.0;
         printf("[%d]: %e %e %e %e\n", i, tc_out[i], ck_out[i], std::fabs(tc_out[i]-ck_out[i])/refval, refval);
      }
   }
#else
   if (enableTChem)
   {
      fprintf(stderr,"TChem requested but not available. Recompile with TChem enabled\n");
      return 2;
   }
#endif

#ifdef USE_PYJAC
   CV::Integrator< PyJac::Functor > cv_pyjac( neq );
   PyJac::Functor pyjac_func;
   if (enablePyJac)
   {
      VectorType<double> pj_out(neq), ck_out(neq);

      cklib_func( neq, 0.0, &u_in[0], ck_out.getPointer());
      pyjac_func( neq, 0.0, &u_in[0], pj_out.getPointer());
      for (int i = 0; i < neq; ++i)
      {
         double refval = std::fabs(pj_out[i]);
         if ( refval <= 1e-20 ) refval = 1.0;
         printf("[%d]: %e %e %e %e\n", i, pj_out[i], ck_out[i], std::fabs(pj_out[i]-ck_out[i])/refval, refval);
      }
# ifdef USE_TCHEM
      if (enableTChem)
      {
         fprintf(stderr,"Error! Both TChem and PyJac were requested. Deselect one to avoid confusion.\n");
         return 10;
      }
# endif
   }
#else
   if (enablePyJac)
   {
      fprintf(stderr,"PyJac requested but not available. Recompile with PyJac enabled\n");
      return 2;
   }
#endif

   VectorType<double> _rwk;
   VectorType<int> _iwk;

   const int stride = 1;
   int nst = 0, nit = 0, nfe = 0, nje = 0;
   double ysum = 0;

   const double time_start = WallClock();

   for (int problem_id = 0; problem_id < numProblems; problem_id += stride)
   {
      double t0 = WallClock();

      for (int k = 0; k < neq; ++k)
         u[k] = u_in[problem_id * neq + k];

      int _nst = 0, _nit = 0, _nfe = 0, _nje = 0;

      double Temp0 = u[ getTempIndex(neq) ];

      if (solverTag == CV)
      {
#ifdef USE_SUNDIALS
         auto _ret = 
#ifdef USE_TCHEM
            ( enableTChem ) ?
               cv_driver (neq, u.getPointer(), t_stop, cv_tchem, tchem_func ) :
#endif
#ifdef USE_PYJAC
            ( enablePyJac ) ?
               cv_driver (neq, u.getPointer(), t_stop, cv_pyjac, pyjac_func ) :
#endif
               cv_driver (neq, u.getPointer(), t_stop, cv_obj, cklib_func );

         _nst += std::get<0>(_ret);
         _nfe += std::get<1>(_ret);
         _nje += std::get<2>(_ret);
         _nit += std::get<3>(_ret);
#else
         std::cerr << "recompile with -DUSE_SUNDIALS" << std::endl;
#endif
      }
      else if (solverTag == RK)
      {
         rk_t rk_;
         rk_counters_t counters_;

         double t_ = 0, h_ = 0;
         rk_create (&rk_, neq);

         rk_.max_iters = 1000;
         rk_.min_iters = 1;

         int _lrwk = rk_lenrwk (&rk_);
         if (_rwk.size() != _lrwk) _rwk.resize(_lrwk);

         rk_init (&rk_, t_, t_stop);

         int ierr_ = rk_solve (&rk_, &t_, &h_, &counters_, &u[0], &_rwk[0], rhs_func, &cklib_func);
         if (ierr_ != RK_SUCCESS)
         {
            fprintf(stderr,"%d: rk_solve error %d %d %d\n", problem_id, ierr_, counters_.niters, rk_.max_iters);
            exit(-1);
         }

         _nst = counters_.nsteps;
         _nit = counters_.niters;
         _nfe = _nit * 6;

         rk_destroy(&rk_);
      }
      else if (solverTag == ROS)
      {
         ros_t ros_;
         ros_counters_t counters_;
         double t_ = 0, h_ = 0;
         ros_create (&ros_, neq, Ros4);

         int _lrwk = ros_lenrwk (&ros_);
         int _liwk = ros_leniwk (&ros_);
         if (_rwk.size() != _lrwk) _rwk.resize(_lrwk);
         if (_iwk.size() != _liwk) _iwk.resize(_liwk);

         ros_init (&ros_, t_, t_stop);

         ros_solve (&ros_, &t_, &h_, &counters_, &u[0], &_iwk[0], &_rwk[0], rhs_func, NULL, &cklib_func);

         _nst = counters_.nst;
         _nit = counters_.niters;
         _nfe = counters_.nfe;
         _nje = counters_.nje;

         ros_destroy(&ros_);
      }
      else if (solverTag == SDIRK)
      {
         sdirk_t sdirk_obj;
         sdirk_counters_t counters_;
         double t_ = 0, h_ = 0;
         sdirk_create (&sdirk_obj, neq, S4a);

         int _lrwk = sdirk_lenrwk (&sdirk_obj);
         int _liwk = sdirk_leniwk (&sdirk_obj);
         if (_rwk.size() != _lrwk) _rwk.resize(_lrwk);
         if (_iwk.size() != _liwk) _iwk.resize(_liwk);

         sdirk_init (&sdirk_obj, t_, t_stop);

         sdirk_solve (&sdirk_obj, &t_, &h_, &counters_, &u[0], &_iwk[0], &_rwk[0], rhs_func, NULL, &cklib_func);

         _nst = counters_.nst;
         _nit = counters_.niters;
         _nfe = counters_.nfe;
         _nje = counters_.nje;

         sdirk_destroy(&sdirk_obj);
      }

      nst += _nst;
      nit += _nit;
      nfe += _nfe;
      nje += _nje;

      double t1 = WallClock();
      int write_stride = std::max(1,numProblems / 16);
      if (problem_id % write_stride == 0)
         printf("%d: %d %d %d %d %f %f %f\n", problem_id, _nst, _nit, _nfe, _nje, u[getTempIndex(neq)], Temp0, 1000*(t1-t0));

      ysum += u[getTempIndex(neq)];
   }

   printf("ysum = %f\n", ysum);
   double calc_time = WallClock() - time_start;
   printf("time = %f %f %d\n", calc_time, calc_time / numProblems, numProblems);
   printf("nst = %d, nit = %d, nfe = %d, nje = %d\n", nst, nit, nfe, nje);

   if (ckptr)
      ck_destroy(&ckptr);

   return 0;
}
