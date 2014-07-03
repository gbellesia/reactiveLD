#include <unordered_map>
#include <vector>
#include <set>

#include "util.hpp"
#include "Reactions.hpp"

class Particle {
public:
  double x, y, z;
  int type;

  //Need this constructor to use this class in maps (with the [] operator)o
  Particle() {};
  Particle(double x, double y, double z, int type) : x(x), y(y), z(z), type(type) {};
};

class Particles {
public:
  std::default_random_engine generator;

  int Nx, Ny, Nz; //Input
  double X, Y, Z;

  double dX, dY, dZ;

  bool spherical, periodic;

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

  Particles(double maxR, double X, double Y, double Z, bool spherical, bool periodic, Reactions &reacs) : X(X), Y(Y), Z(Z), spherical(spherical), periodic(periodic), reacs(reacs)
  {
    Nx = int(X / maxR) + 1;
    Ny = int(Y / maxR) + 1;
    Nz = int(Z / maxR) + 1;

    dX = X / Nx;
    dY = Y / Ny;
    dZ = Z / Nz;
  }

  // Insert a particle of given type
  //   return index if success
  //   print error if failure and quit //return -1 if failure
  int insertParticle(double x, double y, double z, int type, int id = -1)
  {
    indexListT otherParticles = collide(x, y, z, type, id);

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

    particles[id] = Particle(x, y, z, type);

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
    double xx, yy, zz;

    double radius = reacs.bam[type].radius;

    xx = std::abs(adjust(x, X) - X / 2.0) + radius;
    yy = std::abs(adjust(y, Y) - Y / 2.0) + radius;
    zz = std::abs(adjust(z, Z) - Z / 2.0) + radius;

    if(!periodic)
      {
        if((x - radius) < 0 || (x + radius) > X || (y - radius) < 0 || (y + radius) > Y || (z - radius) < 0 || (z + radius) > Z)
          {
            return indexListT({ -1 });
          }
      }
    
    if(spherical)
      {
        if((4 * xx * xx / (X * X) + 4 * yy * yy / (Y * Y) + 4 * zz * zz / (Z * Z)) >= 1.0)
          {
            return indexListT({ -1 });
          }
      }

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

                    if(periodic)
                      {
                        dx = cyclicDistance(pj.x, x, X);
                        dy = cyclicDistance(pj.y, y, Y);
                        dz = cyclicDistance(pj.z, z, Z);
                        }
                    else
                      {
                        dx = std::abs(pj.x - x);
                        dy = std::abs(pj.y - y);
                        dz = std::abs(pj.z - z);
                      }

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
    insertParticle(p.x, p.y, p.z, p.type, id);
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
