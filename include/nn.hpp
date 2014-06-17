#include <unordered_map>
#include <vector>
#include <set>

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

  Particles(int Nx, int Ny, int Nz, double X, double Y, double Z, Reactions &reacs) : Nx(Nx), Ny(Ny), Nz(Nz), X(X), Y(Y), Z(Z), reacs(reacs)
  {
    dX = X / Nx;
    dY = Y / Ny;
    dZ = Z / Nz;
  }

  // Insert a particle of given type
  //   return index if success
  //   print error if failure and quit //return -1 if failure
  int insertParticle(double x, double y, double z, int type, int id = -1)
  {
    indexListT otherParticles = collide(x, y, z, type);

    if(otherParticles.size() > 0)
      {
        std::cout << "Failed to insert new particle. Check collisions first" << std::endl;

        exit(-1);
      }

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

    bins[idx(x, y, z)].insert(id);

    return id;
  }

  // Return a list of all particles that collide with this atom 
  indexListT collide(double x, double y, double z, int type)
  {
    // We use a 3x3x3 collision grid
    // This requires each bin to be 1x times the size of the largest particle

    indexListT indexList;

    std::vector<int> ox = { -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    std::vector<int> oy = { -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1 };
    std::vector<int> oz = { -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1 };

    for(int i = 0; i < (int)ox.size(); i++)
      {
        int idxx = idx(x + ox[i], y + oy[i], z + oz[i]);

        auto binit = bins.find(idxx);
        if(binit != bins.end())
          {
            for(auto it = binit->second.begin(); it != binit->second.end(); it++)
              {
                Particle pj = particles[*it];

                double rsq = (pj.x - x) * (pj.x - x) + (pj.y - y) * (pj.y - y) + (pj.z - z) * (pj.z - z);
  
                if(rsq < (reacs.bam[type].radius + reacs.bam[pj.type].radius) * (reacs.bam[type].radius + reacs.bam[pj.type].radius))
                  {
                    indexList.insert(indexList.end(), binit->second.begin(), binit->second.end());
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
    
    deleteParticle(it->first);

    p.x = x;
    p.y = y;
    p.z = z;

    insertParticle(p.x, p.y, p.z, p.type, id);
  }

  // Removes particle from particle list and bins
  // If fails to delete, just exit program with error
  void deleteParticle(int id)
  {
    auto it = particles.find(id);

    if(it == particles.end())
      {
        std::cout << "Trying to delete a non-existant particle :" << id << std::endl;

        exit(-1);
      }

    Particle &p = it->second;

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

    binIt->second.erase(indexBinIt);
  }
};
