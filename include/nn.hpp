#include <unordered_map>
#include <vector>
#include <set>

#include "util.hpp"
#include "Reactions.hpp"

enum BoxType {
  ellipsoid
};

class Particle {
public:
  double x, y, z;
  double fx, fy, fz;
  double px, py, pz;
  int type;

  //Need this constructor to use this class in maps (with the [] operator)o
  Particle() {};
  Particle(double x, double y, double z, double px, double py, double pz, int type) : x(x), y(y), z(z), px(px), py(py), pz(pz), type(type) {};

  std::string str()
  {
    std::stringstream ss;

    ss << "x : " << x << std::endl <<
      "y : " << y << std::endl <<
      "z : " << z << std::endl <<
      "px : " << px << std::endl <<
      "py : " << py << std::endl <<
      "pz : " << pz << std::endl <<
      "fx : " << fx << std::endl <<
      "fy : " << fy << std::endl <<
      "fz : " << fz << std::endl <<
      "type : " << type << std::endl;

    return ss.str();
  }
};

class Particles {
public:
  std::default_random_engine generator;

  int Nx, Ny, Nz; //Input
  double X, Y, Z;

  double dX, dY, dZ;

  double ll;

  BoxType boxType;

  Reactions &reacs;

  typedef std::vector<int> indexListT;
  typedef std::set<int> indexSetT;
  typedef std::unordered_map< int, indexSetT > binsT;

  binsT bins;
 
  std::unordered_map<int, Particle> particles;

  int idx(double x, double y, double z)
  {
    return idx(int(x / dX), int(y / dY), int(z / dZ));
  }

  int idx(int x, int y, int z)
  {
    return (x % Nx) * Ny * Nz + (y % Ny) * Nz + (z % Nz);
  }

  Particles(double maxR, double X, double Y, double Z, double ll, BoxType boxType, Reactions &reacs) : X(X), Y(Y), Z(Z), ll(ll), boxType(boxType), reacs(reacs)
  {
    Nx = int(X / (3.0 * maxR)) + 1;
    Ny = int(Y / (3.0 * maxR)) + 1;
    Nz = int(Z / (3.0 * maxR)) + 1;

    dX = X / Nx;
    dY = Y / Ny;
    dZ = Z / Nz;
  }

  // Insert a particle of given type
  //   return index if success
  //   print error if failure and quit //return -1 if failure
  int insertParticle(double x, double y, double z, double px, double py, double pz, int type, int id = -1)
  {
    //indexListT otherParticles = collide(x, y, z, type, id);

    /*if(otherParticles.size() > 0)
      {
        std::cout << "Failed to insert new particle. Check collisions first" << std::endl;

        exit(-1);
        }*/

    if(id == -1)
      {
        id = generator();
        while(particles.find(id) != particles.end())
          {
            id = generator();
          }
      }
    else
      {
        if(particles.find(id) != particles.end())
          {
            std::cout << "Requested id " << id << " already exists in particles list" << std::endl;
          }
      }

    particles[id] = Particle(x, y, z, px, py, pz, type);

    //std::cout << x << " " << y << " " << z << " " << std::endl;
    //std::cout << X << " " << Y << " " << Z << " " << std::endl;

    int bin = idx(x, y, z);

    //std::cout << "Inserting " << id << " at " << bin << std::endl;

    bins[bin].insert(id);

    auto it = bins[bin].find(id);

    if(it == bins[bin].end())
      std::cout << "ERRORORORORORORORORORORORORRRRRRRRRR" << std::endl;

    return id;
  }

  // Return a list of all particles that collide with this atom 
  indexListT collide(double x, double y, double z, int type, int id = -1)
  {
    // We use a 3x3x3 collision grid
    // This requires each bin to be 1x times the size of the largest particle

    indexListT indexList;

    std::vector<int> ox = { -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    std::vector<int> oy = { -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1 };
    std::vector<int> oz = { -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1 };
    //std::vector<int> ox = { 0 };
    //std::vector<int> oy = { 0 };
    //std::vector<int> oz = { 0 };

    std::map<int, int> test;
    std::vector<int> visted;

    for(int i = 0; i < (int)ox.size(); i++)
      {
        int idxx = idx(int(x / dX) + ox[i], int(y / dY) + oy[i], int(z / dZ) + oz[i]);
        visted.push_back(idxx);

        auto binit = bins.find(idxx);
        if(binit != bins.end())
          {
            for(auto it = binit->second.begin(); it != binit->second.end(); it++)
              {
                if(id != *it)
                {
                  if(test.find(*it) != test.end())
                    {
                      std::cout << "panic" << std::endl;
                    }
                  
                  Particle pj = particles[*it];

                  double dx, dy, dz;
                  
                  dx = std::abs(pj.x - x);
                  dy = std::abs(pj.y - y);
                  dz = std::abs(pj.z - z);
                  
                  double rsq = dx * dx + dy * dy + dz * dz;
                  
                  if(rsq < (reacs.bam[type].radius + reacs.bam[pj.type].radius) * (reacs.bam[type].radius + reacs.bam[pj.type].radius))
                    {
                      indexList.push_back(*it);
                      test[*it] = idxx;
                      //indexList.insert(indexList.end(), binit->second.begin(), binit->second.end());
                    }
                }
              }
          }
      }

    return indexList;
  }

  void computeForces(Reactions &reacs)
  {
    // Initialize forces to boundary forces
    for(auto it = particles.begin(); it != particles.end(); it++)
      {
        Particle &particle = it->second;

        double k1 = 1.0;
        double k2 = -1.0;
        double exp1 = 2.0; // exponent of the potential
        double exp2 = 2.0; // exponent of the potential
        double dist2 = (particle.x - X / 2.0) * (particle.x - X / 2.0) + (particle.y - Y / 2.0) * (particle.y - Y / 2.0) + (particle.z - Z / 2.0) * (particle.z - Z / 2.0);
        double dist = sqrt(dist2);
        double R = reacs.bam[particle.type].radius;
        double L = X / 2.0;
        double cx, cy, cz;

        particle.fx = 0;
        particle.fy = 0;
        particle.fz = 0;

        double fval, rval;
        if (dist > L) {
          cx = (particle.x - X / 2.0) / dist;
          cy = (particle.y - Y / 2.0) / dist;
          cz = (particle.z - Z / 2.0) / dist;
          rval = fabs(dist - L);
          fval = - (exp1) * k1;
          fval *= pow(rval, exp1 - 1);
          particle.fx = fval * cx;
          particle.fy = fval * cy;
          particle.fz = fval * cz;
        } else {
          particle.fx = 0.0;
          particle.fy = 0.0;
          particle.fz = 0.0;
        }

        /*if (dist > (L - R)) {
          cx = (particle.x - X / 2.0) / dist;
          cy = (particle.y - Y / 2.0) / dist;
          cz = (particle.z - Z / 2.0) / dist;
          rval = fabs(dist - (L-R));
          fval = - (exp2) * k2;
          fval *= pow(rval, exp2-1);
          particle.fx -= fval * cx;
          particle.fy -= fval * cy;
          particle.fz -= fval * cz;
        } else {
          particle.fx += 0.0;
          particle.fy += 0.0;
          particle.fz += 0.0;
          }*/

        /*std::cout << particle.str() << std::endl;
        std::cout << "dist " << dist << std::endl;
        std::cout << std::endl;*/
      }

    for(auto binit = bins.begin(); binit != bins.end(); binit++)
      {
        for(auto it = binit->second.begin(); it != binit->second.end(); it++)
          {
            for(auto it2 = std::next(it); it2 != binit->second.end(); it2++)
              {
                Particle &particle1 = particles[*it], &particle2 = particles[*it2];

                // This adds the incremental bit of force between these two particles
                double dx = particle1.x - particle2.x,
                  dy = particle1.y - particle2.y,
                  dz = particle1.z - particle2.z;

                double r2 = sqrt(dx * dx + dy * dy + dz * dz);
                double sigma = (reacs.bam[particle1.type].radius + reacs.bam[particle2.type].radius) / 2.0;
                double eps = sqrt(reacs.bam[particle1.type].eps * reacs.bam[particle2.type].eps);
                double rc2 = pow(3.0*sigma,2);

                double ljrcut = 24.0 * eps *(2.0 * pow(sigma,12)/pow(rc2,7) - ll * pow(sigma,6)/pow(rc2,4));

                double f = 24.0 * eps *(2.0 * pow(sigma,12)/pow(r2,7) - ll * pow(sigma,6)/pow(r2,4)) - ljrcut;

                particle1.fx += f * dx;
                particle1.fy += f * dy;
                particle1.fz += f * dz;

                particle2.fx -= f * dx;
                particle2.fy -= f * dy;
                particle2.fz -= f * dz;
              }
          }
      }
  }

  // Indicate that particle has moved. Update it in the index. Keep the same id
  void move(int id, double x, double y, double z)
  {
    auto it = particles.find(id);

    if(it == particles.end())
      {
        std::cout << "Trying to move a non-existant particle :" << id << std::endl;

        exit(-1);
      }

    Particle p = it->second;
    
    //printf("Movin', yo!");
    deleteParticle(it->first);

    p.x = x;
    p.y = y;
    p.z = z;

    //std::cout << "inserting particle " << 
    insertParticle(p.x, p.y, p.z, p.px, p.py, p.pz, p.type, id);
    // << std::endl;
  }

  // Removes particle from particle list and bins
  // If fails to delete, just exit program with error
  void deleteParticle(int id)
  {
    auto it = particles.find(id);

    //std::cout << "Deleting particle" << id << "\n";

    if(it == particles.end())
      {
        std::cout << "Trying to delete a non-existant particle :" << id << std::endl;

        exit(-1);
      }

    Particle &p = it->second;

    //std::cout << " from bin " << idx(p.x, p.y, p.z) << std::endl;

    auto binIt = bins.find(idx(p.x, p.y, p.z));

    if(binIt == bins.end())
      {
        std::cout << "Particle is not in the right bin (when deleting)\n";
        
        exit(-1);
      }

    auto indexBinIt = binIt->second.find(id);

    if(indexBinIt == binIt->second.end())
      {
        std::cout << "Particle was not in expected bin (when del.)\n";

        exit(-1);
      }

    particles.erase(it);

    binIt->second.erase(indexBinIt);
  }
};
