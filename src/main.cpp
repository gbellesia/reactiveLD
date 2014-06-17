#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <random>
#include "functions.h"
#include <string.h>
#include <sstream>
#include <random>
#include <nn.hpp>
#include <algorithm>
#include <set>
#include <jansson.h>

#include "Reactions.hpp"

#define debug false

int main() {
  std::random_device rd;
  std::mt19937 generator(rd());
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::normal_distribution<double> gaussian(0.0, 1.0);
  std::stringstream ss;

  double X = 10, Y = 10, Z = 10;
  double D = 0.5;

  int N = 10;

  int maxTries = 1000;

  Reactions reacs;

  double dt = 1.0;
  double steps = 100;

  // Read in configuration file
  std::string lineInput;
  while (std::getline(std::cin, lineInput)) {
    ss << lineInput;
  }

  json_error_t error;

  json_t *root = json_loads(ss.str().c_str(), 0, &error);

  if(!root)
    {
      fprintf(stderr, "error: on line %d: %s\n", error.line, error.text);
      return 1;
    }

  if(!json_is_object(root))
    {
      fprintf(stderr, "error: root is not an object\n");
      json_decref(root);
      return 1;
    }

  json_t *mArray = json_object_get(root, "monatomic");
  json_t *bArray = json_object_get(root, "binary");
  json_t *aTypes = json_object_get(root, "types");

  for(int i = 0; i < (int)json_array_size(bArray); i++)
    {
      json_t *r = json_array_get(bArray, i);
      if(!json_is_object(r))
        {
          fprintf(stderr, "Found a non-object item where binary reactions should be\n");
          json_decref(root);
          return 1;
        }

      json_t *Aj = json_object_get(r, "A");
      if(Aj == NULL || !json_is_integer(Aj))
        {
          std::cout << "Found a non-number identifier for (A) atoms" << std::endl;

          if(Aj != NULL)
            json_decref(root);

          return 1;
        }

      json_t *Bj = json_object_get(r, "B");
      if(Bj == NULL || !json_is_integer(Bj))
        {
          std::cout << "Found a non-number identifier for B-atoms" << std::endl;

          if(Bj != NULL)
            json_decref(root);

          return 1;
        }

      json_t *Cj = json_object_get(r, "C");
      if(Cj == NULL || !json_is_integer(Cj))
        {
          std::cout << "Found a non-number identifier for (C) atoms" << std::endl;

          if(Cj != NULL)
            json_decref(root);

          return 1;
        }

      json_t *kj = json_object_get(r, "k");
      if(kj == NULL || !json_is_real(kj))
        {
          std::cout << "Found a non-real identifier for reaction rates" << std::endl;

          if(kj != NULL)
            json_decref(root);

          return 1;
        }

      int C = json_integer_value(Cj);
      int B = (Bj == NULL) ? -1 : json_integer_value(Bj);
      int A = json_integer_value(Aj);
     
      double k = json_real_value(kj);

      reacs.brm.insert(Reactions::brmT::value_type(mpp(A, B), BinaryReaction(A, B, C, k)));
    }

  for(int i = 0; i < (int)json_array_size(mArray); i++)
    {
      json_t *r = json_array_get(mArray, i);
      if(!json_is_object(r))
        {
          fprintf(stderr, "Found a non-object item where monatomic reactions should be\n");
          json_decref(root);
          return 1;
        }

      json_t *Aj = json_object_get(r, "A");
      if(Aj == NULL || !json_is_integer(Aj))
        {
          std::cout << "Found a non-number identifier for (A) atoms" << std::endl;

          if(Aj != NULL)
            json_decref(root);

          return 1;
        }

      json_t *Bj = json_object_get(r, "B");
      if(Bj == NULL || !json_is_integer(Bj))
        {
          std::cout << "Found a non-number identifier for B-atoms" << std::endl;

          json_decref(root);

          return 1;
        }

      json_t *Cj = json_object_get(r, "C");
      if(Cj != NULL && !json_is_integer(Cj))
        {
          std::cout << "Found a non-number identifier for (C) atoms" << std::endl;

          if(Cj != NULL)
            json_decref(root);

          return 1;
        }

      json_t *kj = json_object_get(r, "k");
      if(kj == NULL || !json_is_real(kj))
        {
          std::cout << "Found a non-real identifier for reaction rates" << std::endl;

          if(kj != NULL)
            json_decref(root);

          return 1;
        }

      int C = json_integer_value(Cj);
      int B = (Bj == NULL) ? -1 : json_integer_value(Bj);
      int A = json_integer_value(Aj);
     
      double k = json_real_value(kj);

      reacs.mrm.insert(Reactions::mrmT::value_type(A, MonatomicReaction(A, B, C, k)));
    }

  for(int i = 0; i < (int)json_array_size(aTypes); i++)
    {
      json_t *r = json_array_get(aTypes, i);
      if(!json_is_object(r))
        {
          fprintf(stderr, "Found a non-object item where monatomic reactions should be\n");
          json_decref(root);
          return 1;
        }

      json_t *typeJ = json_object_get(r, "type");
      if(typeJ == NULL || !json_is_integer(typeJ))
        {
          std::cout << "Found a non-number identifier for (A) atoms" << std::endl;

          if(typeJ != NULL)
            json_decref(typeJ);

          return 1;
        }

      json_t *radiusJ = json_object_get(r, "radius");
      if(radiusJ == NULL || !json_is_real(radiusJ))
        {
          std::cout << "Found a non-real identifier for reaction rates" << std::endl;

          if(radiusJ != NULL)
            json_decref(radiusJ);

          return 1;
        }

      json_t *Dj = json_object_get(r, "D");

      if(Dj == NULL || !json_is_real(Dj))
        {
          std::cout << "Found a non-number identifier for B-atoms" << std::endl;

          if(Dj != NULL)
            json_decref(Dj);

          return 1;
        }

      int type = json_integer_value(typeJ);
      double radius = json_real_value(radiusJ);
      double D = json_real_value(Dj);

      reacs.bam[type] = BrownianAtom(type, radius, D);
    }

  reacs.init(dt);

  if(reacs.bam.size() < 1)
    {
      std::cout << "Not enough atoms to make a simulation!" << std::endl;

      return 1;
    }

  Particles parts(10, 10, 10, X, Y, Z, reacs);

  // Initialize particles
  for(int n = 0; n < N; n++)
    {
      for(int tt = 0; tt < maxTries; tt++)
        {
          double x = X * uniform(generator),
            y = Y * uniform(generator),
            z = Z * uniform(generator);

          Particles::indexListT indexList = parts.collide(x, y, z, 1);

          if(indexList.size() < 1)
            {
              parts.insertParticle(x, y, z, 1);

              break;
            }
        }
    }

  for(int s = 0; s < steps; s++)
    {
      std::vector<int> oldParticles;
      std::set<int> deleted;
      
      for(auto it = parts.particles.begin(); it != parts.particles.end(); it++)
        {
          oldParticles.push_back(it->first);
        }

      for(auto it = oldParticles.begin(); it != oldParticles.end(); it++)
        {
          // Don't try to work with an atom that has already been deleted
          if(deleted.find(*it) != deleted.end())
            continue;

          int pid = *it;

          Particle &p = parts.particles[pid];

          double dx = gaussian(generator) * sqrt(2 * D * dt);
          double dy = gaussian(generator) * sqrt(2 * D * dt);
          double dz = gaussian(generator) * sqrt(2 * D * dt);

          double nx = p.x + dx, ny = p.y + dy, nz = p.z + dz;

          //vector of ints pointing to other particles
          auto touching = parts.collide(p.x + dx, p.y + dy, p.z + dz, p.type);

          // Reject the situation where three particles meet
          if(touching.size() > 2 || D == 0.0)
            {
              continue;
            }
          else if(touching.size() == 0) // If there are no collisions, check if this particle can react in any way
            // Any way being: A -> B
            // or   A -> B + C
            {
              std::pair<Reactions::mrmT::iterator, Reactions::mrmT::iterator> se = reacs.mrm.equal_range(p.type);

              Reactions::mrmT::iterator it = select<Reactions::mrmT>(se.first, se.second, reacs.mrm.end(), dt, 1.0, uniform(generator), uniform(generator));

              // If no reactions are found for this particle, accept the move
              if(it == reacs.mrm.end())
                {
                  parts.move(it->first, p.x + dx, p.y + dy, p.z + dz);

                  continue;
                }

              MonatomicReaction reaction = it->second;

              //std::cout << reaction.str() << std::endl;

              int Btype = reaction.B,
                Ctype = reaction.C;

              //If there is no Btype, better be no Ctype. This is a destruction reac!
              if(Btype > 0)
                {
                  int tries = 0;

                  double xx, yy, zz;

                  for(tries = 0; tries < maxTries; tries++)
                    {
                      //Check if it will be possible to insert a B atom
                      auto touching = parts.collide(nx, ny, nz, Btype);

                      // The newly inserted atom might hit the two atoms we're about to delete
                      //   But if it hits a third, or it hits an atom other than the one we're possibly
                      //     removing, reject
                      std::cout << touching.size() << std::endl;
                      if(touching.size() >= 1)
                        {
                          if(pid != touching[0] || touching.size() > 1)
                            {
                              if(debug)
                                printf("Failed to insert new atom, atom still moved \n");
              
                              continue;
                            }
                        }
          
                      //Check if it will be possible to insert a C atom
                      if(Ctype > 0)
                        {
                          double r = reacs.pseps[mpp(Btype, Ctype)].sample(uniform(generator));

                          //std::cout << r << std::endl;
                  
                          double theta = pi * uniform(generator);
                          double phi = 2 * pi * uniform(generator);
                  
                          xx = r * sin(theta) * cos(phi);
                          yy = r * sin(theta) * sin(phi);
                          zz = r * cos(theta);

                          //printf("%f %f %f -> %f %f %f, %f\n", x[atomi][0], x[atomi][1], x[atomi][2], xx, yy, zz, r);
                          auto touching = parts.collide(p.x + xx, p.y + yy, p.z + zz, Ctype);

                          if(touching.size() >= 1)
                            {
                              if(pid != touching[0] || touching.size() > 1)
                                {
                                  if(debug)
                                    printf("Failed to insert new atom, atom still moved \n");
                      
                                  continue;
                                }
                            }
                        }

                      // If we got this far, any necessary reaction products have been successfully inserted
                      break;
                    }

                  // If we failed to insert the new atoms, reject the move and go on
                  if(tries == maxTries)
                    {
                      //std::cout << "reject" << std::endl;
              
                      continue;
                    }
                  else
                    {                      
                      // Get rid of the reacted particle
                      deleted.insert(pid);

                      parts.deleteParticle(pid);
                  
                      // Insert the two products
                      parts.insertParticle(nx, ny, nz, Btype);

                      if(Ctype > 0)
                        {
                          parts.insertParticle(nx + xx, ny + yy, nz + zz, Ctype);
                        }
                    }
                }
            }
          else // Handle binary reactions
            {
              // Remember that findTouching call way above? Let's get the atom id out of that
              int pjd = touching[0];

              Particle &p2 = parts.particles[pjd];

              // We have an i touching a j, so build the appropriate particle pair
              ParticlePair pp = mpp(p.type, parts.particles[pjd].type);
          
              // We could either have a regular binary reaction or an enzymatic reaction
              std::pair<Reactions::brmT::iterator, Reactions::brmT::iterator> se = reacs.brm.equal_range(pp);

              Reactions::brmT::iterator it = select<Reactions::brmT>(se.first, se.second, reacs.brm.end(), dt, reacs.paccs[pp], uniform(generator), uniform(generator));

              //if(type[atomi] != type[pjd])
              //  printf("particle %d and %d have collided!\n", type[atomi], type[pjd]);

              // If no reactions are found for this particle, reject the move, and end
              if(it == reacs.brm.end())
                {
                  //if(type[atomi] != type[pjd])
                  //printf("No reaction found!\n");
              
                  continue;
                }
              //printf("Reaction\n");

              BinaryReaction reaction = it->second;

              double nx = p2.x, ny = p2.y, nz = p2.z;

              int Ctype = reaction.C;

              if(Ctype > 0)
                {
                  //Check if it will be possible to insert a C atom
                  auto touching = parts.collide(nx, ny, nz, Ctype);

                  std::set<int> tmp = { pid, pjd };
                  if(touching.size() > 2 ||
                     !std::includes(tmp.begin(), tmp.end(), touching.begin(), touching.end()))
                    {
                      if(debug)
                        printf("Failed to insert new atom \n");
              
                      continue;
                    }
                }
          
              // Time to delete the reactants
              if(debug)
                printf("Trying to delete\n");

              deleted.insert(pid);
              deleted.insert(pjd);

              parts.deleteParticle(pid);
              parts.deleteParticle(pjd);

              parts.insertParticle(nx, ny, nz, Ctype);
            }
        }
    }
}
