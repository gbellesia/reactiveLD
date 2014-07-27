#ifdef __REACTIONS_HPP__
#else

#define __REACTIONS_HPP__

#include <set>
#include <map>
#include <vector>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <string>
#include <sstream>

const double pi = 3.141592653589793238462643383279502884;

template<class T>
typename T::iterator select(typename T::iterator it1, typename T::iterator it2, typename T::iterator itend, double dt, double pacc, std::uniform_real_distribution<double> &uniform, std::mt19937 &generator)
{
  std::vector<typename T::iterator> opts;
  
  for(typename T::iterator it = it1; it != it2; it++)
    {
      // Check propensity of firing:
      double r = uniform(generator);
      bool fire = r < it->second.k * dt * pacc;

      //printf("A %d B %d r = %f, %f, %f\n", it->second.A, it->second.B, r, it->second.k * dt * pacc, pacc);

      // If the reaction doesn't fire, accept the move and we're done
      if(!fire)
        {
          continue;
        }
      
      opts.push_back(it);

      //std::cout << i++ << std::endl;
    }

  double sum = 0.0;

  if(opts.size() == 0)
    return itend;

  for(int i = 0; i < (int)opts.size(); i++)
    sum += (*opts[i]).second.k * dt * pacc;
  
  double v = uniform(generator) * sum;
  double tsum = 0.0;

  //std::cout << "v" << v << " " << sum << std::endl;
  for(int i = 0; i < (int)opts.size(); i++)
    {
      tsum += (*opts[i]).second.k * dt * pacc;
      if(tsum >= v)
	{
	  return opts[i];
	}
    }
  
  return itend;
}

// This is the CDF that is used to sample the inverse radial distributions
class PsepCDF {
public:
  struct FuncParams
  {
    double R, b, N;
  };
  
  static double func(double x, void *params)
  {
    FuncParams &p = *((FuncParams *)params);
    
    return (4 * pi * p.R * x / p.N) * erfc((x - p.R) / p.b);
  };

  const static int bins = 1000;

  double R, D, dt;
  std::map<double, double> steps;

  FuncParams fp;

  //I think this needs to be here for the maps and such to work
  PsepCDF() {};
  PsepCDF(double R, double D, double dt) : R(R), D(D), dt(dt) {
    fp.b = sqrt(4 * D * dt);
    fp.N = 4 * pi * (R * D * dt + 2 * R * R * sqrt(D * dt / pi));
    fp.R = R;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

    double result, error;

    gsl_function F;
    F.function = &func;
    F.params = (void *)&fp;
    //printf("%f\n", func(fp.R, &fp));
    double dy = 20 * fp.R / bins;
    double lasty = -1.0;

    //The random number will be between 0.0 and 1.0, gotta have all these points!
    steps[0.0] = fp.R;
    for(int i = 0; i < bins; i++)
      {
        gsl_integration_qags (&F, i * dy + fp.R, (i + 1) * dy + fp.R, 0, 1e-7, 1000,
                              w, &result, &error);

        if(i != 0)
          {
            result = result + lasty;
          }

        //printf("result %f\n", result);
        
        //Gotta have the 1.0 point as well!
        if(i == bins - 1)
          {
            steps[1.0] = (i + 1) * dy + fp.R;
          }
        else
          {
            steps[result] = (i + 1) * dy + fp.R;
          }

        lasty = result;
      }

    gsl_integration_workspace_free (w);
  }

  double sample(double r)
  {
    std::map<double, double>::iterator itr = steps.lower_bound(r), itl;

    //printf("%f %d %f %f\n", r, steps.size(), steps.begin()->first, steps.begin()->second);

    // This is technically possible but unlikely. We'll probably need to do some interpolation
    if(itr->first == r)
      return itr->second;

    if(itr == steps.begin())
      {
        printf("SOMETHING IS FISHY IN THE CDF SAMPLING\n");
        exit(-1);
      }

    itl = itr;
    itl--;

    //printf("itl %f %f, itr %f %f\n", itl->first, itl->second, itr->first, itr->second);

    /*for(std::map<double, double>::iterator it = steps.begin(); it != steps.end(); it++)
      {
      printf("%f %f\n", it->first, it->second);
      }*/

    return (itr->first - r) / (itr->first - itl->first) * (itr->second - itl->second) * -1 + itr->second;
  }
};

class BrownianAtom {
public:
  int type;
  double radius, D, mass, eps;
    
  BrownianAtom() {};
  BrownianAtom(int type, double radius, double D, double mass, double eps) : type(type), radius(radius), D(D), mass(mass), eps(eps) {};

  std::string str()
  {
    std::stringstream out;

    out << "Atom type: " << type << ", D: " << D << ", R: " << radius << ", mass: " << mass << ", eps: " << eps;

    return out.str();
  }
};

// For the two reaction classes, any optional species not used will be marked with a -1

// A + B > C
// A + E -> B + E
// A + E -> B + C + E // Not implemented and I won't implement it without a reason :D:D:D:D:D
class BinaryReaction {
public:
  int A, B, C, D;
  double k;
  bool enzymatic;

  BinaryReaction() {};
  BinaryReaction(int A, int B, int C, int D, double k, bool enzymatic = false) : A(A), B(B), C(C), D(D), k(k), enzymatic(enzymatic) {};

  std::string str()
  {
    std::stringstream out;

    out << "Binary reac. A: " << A << ", B: " << B << ", C: " << C << ", D: " << D << ", rate: " << k << ", enzymatic: " << (enzymatic ? "true" : "false");

    return out.str();
  }
};
  
// A -> B
// A -> B + C
class MonatomicReaction {
public:
  int A, B, C;
  double k;

  MonatomicReaction() {};
  MonatomicReaction(int A, int B, int C, double k) : A(A), B(B), C(C), k(k) {}

  std::string str()
  {
    std::stringstream out;

    out << "Monatomic reac. A: " << A << ", B: " << B << ", C: " << C << ", rate: " << k;

    return out.str();
  }
};

//Particle Combos
typedef std::set<int> ParticlePair;
  
// This is a helper function -- "Make Particle Pair" -- because we aren't
//   using a full, carefully defined class for our ParticlePair
ParticlePair mpp(int i, int j)
{
  int array[2] = {i, j};

  return ParticlePair(array, array + 2);
};

class Type {
public:
  double radius;
  double D;

  Type(double radius, double D) : radius(radius), D(D) {}
};

class Reactions {
public:
  //template<class T>
  //typename T::iterator select(typename T::iterator it1, typename T::iterator it2, typename T::iterator itend, double pacc);


  //List of all possible Pseprs
  typedef std::map< ParticlePair, PsepCDF > PsepsT;
  PsepsT pseps;
  //std::map< ParticlePair, PcolCDF > Pcols; NOT IMPLEMENTED YET -- PART OF FRAZIER AND ALBER PAPER
  typedef std::map< ParticlePair, double > PaccsT;
  PaccsT paccs;

  typedef std::unordered_map<int, Type> typesT;
  typesT types;

  typedef std::map< int, BrownianAtom > bamT;
  bamT bam; // Brownian atom map - BAM!

  typedef std::multimap< ParticlePair, BinaryReaction > brmT;
  brmT brm; // Binary reaction map - BRM!

  typedef std::multimap< int, MonatomicReaction > mrmT;
  mrmT mrm; // Monatomic reaction map - MRM!

  // Call this after the types and reactions have all been filled out
  double init(double dt)
  {
    double maxR = 0.0;

    std::cout << "Brownian Atom types" << std::endl;
    for(bamT::iterator it = bam.begin(); it != bam.end(); it++)
      {
        std::cout << it->second.str() << std::endl;
      }
    
    std::cout << std::endl << "Binary stuff" << std::endl;
    
    for(brmT::iterator it = brm.begin(); it != brm.end(); it++)
      {
        std::cout << it->second.str() << std::endl;
      }
    
    std::cout << std::endl << "Monatomic stuff" << std::endl;
    
    for(mrmT::iterator it = mrm.begin(); it != mrm.end(); it++)
      {
        std::cout << it->second.str() << std::endl;
      }
    
    std::cout << std::endl;
    
    for(bamT::iterator it = bam.begin(); it != bam.end(); it++)
    {
      bamT::iterator it2 = it;
      for(; it2 != bam.end(); it2++)
        {
          double D = it->second.D + it2->second.D,
            R = it->second.radius + it2->second.radius;

          maxR = std::max(R, maxR);

          if(D < 1e-17)
            {
              printf("WARNING D < 1e-17 in Psepr evaluations of particle types %d and %d. This is possibly okay if they are immobile.\n", it->first, it2->first);
              continue;
            }

          //printf("%f %f %f\n", D, R, update->dt);

          pseps[mpp(it->first, it2->first)] = PsepCDF(R, D, dt);

          //double N = 4 * pi * (R * D * dt + 2 * R * R * sqrt(D * dt / pi));

          paccs[mpp(it->first, it2->first)] = 1 / pseps[mpp(it->first, it2->first)].fp.N;
        }

      paccs[mpp(it->first, -1)] = 1.0;
    }
    
    return maxR;
  }
};

#endif
