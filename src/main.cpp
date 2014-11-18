#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <random>
#include <string.h>
#include <sstream>
#include <random>
#include <nn.hpp>
#include <algorithm>
#include <set>
#include <jansson.h>
#include <iomanip>

#include "util.hpp"
#include "Reactions.hpp"

#define debug false

int main(int argc, char **argv) {
  double X = 10000, Y = 10000, Z = 10000;

  //int N = 10;

  int maxTries = 1000;

  Reactions reacs;

  double dt;
  int steps;
  int printSteps;

  // Args are [config] [output]

  if(argc < 4)
    {
      std::cout << "Not enough arguments. Input format is: ./bd_run [infile] [outfile] [randomseed]" << std::endl;
    }

  std::string outputType;
  if(argc == 5)
    {
      outputType = std::string(argv[4]);
    }
  else
    {
      outputType = "xyzwhatever";
    }

  std::string inFileName(argv[1]), outFileName(argv[2]);

  std::mt19937 generator(atoi(argv[3]));
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::normal_distribution<double> gaussian(0.0, 1.0);
  std::stringstream ss;

  std::ifstream in(inFileName);

  // Read in configuration file
  std::string lineInput;
  while (std::getline(in, lineInput)) {
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

  json_t *atomsArray = json_object_get(root, "atoms");
  json_t *atomsDistributionArray = json_object_get(root, "atomsDist");
  json_t *mArray = json_object_get(root, "monatomic");
  json_t *Xj = json_object_get(root, "X");
  json_t *Yj = json_object_get(root, "Y");
  json_t *Zj = json_object_get(root, "Z");
  json_t *stepsj = json_object_get(root, "steps");
  json_t *printStepsj = json_object_get(root, "printSteps");
  json_t *dtj = json_object_get(root, "dt");
  json_t *bArray = json_object_get(root, "binary");
  json_t *aTypes = json_object_get(root, "types");
  json_t *type = json_object_get(root, "simType");

  BoxType boxType;
  bool twoD = false;

  if(type == NULL || !json_is_string(type))
    {
      boxType = BoxType::periodicBox;
    }
  else if(strcmp(json_string_value(type), "ellipsoid") == 0)
    {
      boxType = BoxType::ellipsoid;
    }
  else if(strcmp(json_string_value(type), "capsule") == 0)
    {
      boxType = BoxType::capsule;
    }
  else if(strcmp(json_string_value(type), "cylinder") == 0)
    {
      boxType = BoxType::cylinder;
    }
  else if(strcmp(json_string_value(type), "boxWithWalls") == 0)
    {
      boxType = BoxType::boxWithWalls;
    }
  else if(strcmp(json_string_value(type), "periodicBox") == 0)
    {
      boxType = BoxType::periodicBox;
    }
  else
    {
      std::cout << "simType must be defined as either \"ellipsoid\", or \"capsule\" or \"cylinder\", or \"boxWithWalls\", or \"periodicBox\"" << std::endl;
      
      if(Xj != NULL)
        json_decref(root);
      
      return 1;
    }
  
  if(Xj == NULL || !json_is_number(Xj))
    {
      std::cout << "X must be a number" << std::endl;
      
      if(Xj != NULL)
        json_decref(root);
      
      return 1;
    }

  if(Yj == NULL || !json_is_number(Yj))
    {
      std::cout << "Y must be a number" << std::endl;
      
      if(Yj != NULL)
        json_decref(root);
      
      return 1;
    }

  if(Zj == NULL)
    {
      std::cout << "No Z coordinate found. Enforcing 2D mode!" << std::endl;
      twoD = true;
    }

  if(!twoD && !json_is_number(Zj))
    {
      std::cout << "Z must be a number" << std::endl;
      
      if(Zj != NULL)
        json_decref(root);
      
      return 1;
    }

  if(stepsj == NULL || !json_is_integer(stepsj))
    {
      std::cout << "steps must be an integer" << std::endl;
      
      if(stepsj != NULL)
        json_decref(root);
      
      return 1;
    }

  if(printStepsj == NULL || !json_is_integer(printStepsj))
    {
      std::cout << "printSteps must be an integer" << std::endl;
      
      if(printStepsj != NULL)
        json_decref(root);
      
      return 1;
    }
  
  if(dtj == NULL || !json_is_number(dtj))
    {
      std::cout << "dt must be a number" << std::endl;
      
      if(dtj != NULL)
        json_decref(root);
      
      return 1;
    }

  X = json_number_value(Xj);
  Y = json_number_value(Yj);
  Z = (twoD) ? 1.0 : json_number_value(Zj);
  steps = json_integer_value(stepsj);
  printSteps = json_integer_value(printStepsj);
  dt = json_number_value(dtj);

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
      if(Cj != NULL && !json_is_integer(Cj))
        {
          std::cout << "Found a non-number identifier for (C) atoms" << std::endl;

          json_decref(root);

          return 1;
        }

      json_t *Dj = json_object_get(r, "D");
      if(Dj != NULL && !json_is_integer(Dj))
        {
          std::cout << "Found a non-number identifier for (D) atoms" << std::endl;

          json_decref(root);

          return 1;
        }

      json_t *kj = json_object_get(r, "k");
      if(kj == NULL || !json_is_number(kj))
        {
          std::cout << "Found a non-real identifier for reaction rates" << std::endl;

          if(kj != NULL)
            json_decref(root);

          return 1;
        }

      int D = (Dj == NULL) ? -1 : json_integer_value(Dj);
      int C = (Cj == NULL) ? -1 : json_integer_value(Cj);
      int B = json_integer_value(Bj);
      int A = json_integer_value(Aj);
     
      double k = json_number_value(kj);

      reacs.brm.insert(Reactions::brmT::value_type(mpp(A, B), BinaryReaction(A, B, C, D, k)));
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
      if(Aj != NULL && !json_is_integer(Aj))
        {
          std::cout << "Found a non-number identifier for (A) atoms" << std::endl;

          json_decref(root);

          return 1;
        }

      json_t *Bj = json_object_get(r, "B");
      if(Bj != NULL && !json_is_integer(Bj))
        {
          std::cout << "Found a non-number identifier for B-atoms" << std::endl;

          json_decref(root);

          return 1;
        }

      json_t *Cj = json_object_get(r, "C");
      if(Cj != NULL && !json_is_integer(Cj))
        {
          std::cout << "Found a non-number identifier for (C) atoms" << std::endl;

          json_decref(root);

          return 1;
        }

      json_t *kj = json_object_get(r, "k");
      if(kj == NULL || !json_is_number(kj))
        {
          std::cout << "Found a non-real identifier for reaction rates" << std::endl;

          if(kj != NULL)
            json_decref(root);

          return 1;
        }

      int A = (Aj == NULL) ? -1 : json_integer_value(Aj);
      int B = (Bj == NULL) ? -1 : json_integer_value(Bj);
      int C = (Cj == NULL) ? -1 : json_integer_value(Cj);
     
      double k = json_number_value(kj);

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
      if(radiusJ == NULL || !json_is_number(radiusJ))
        {
          std::cout << "Found a non-real identifier for reaction rates" << std::endl;

          if(radiusJ != NULL)
            json_decref(radiusJ);

          return 1;
        }

      json_t *Dj = json_object_get(r, "D");

      if(Dj == NULL || !json_is_number(Dj))
        {
          std::cout << "Found a non-number identifier for B-atoms" << std::endl;

          if(Dj != NULL)
            json_decref(Dj);

          return 1;
        }

      int type = json_integer_value(typeJ);
      double radius = json_number_value(radiusJ);
      double D = json_number_value(Dj);

      reacs.bam[type] = BrownianAtom(type, radius, D);
    }

  double maxR = reacs.init(dt);

  if(reacs.bam.size() < 1)
    {
      std::cout << "Not enough atoms to make a simulation!" << std::endl;

      return 1;
    }

  Particles parts(maxR, X, Y, Z, boxType, twoD, reacs);

  if(atomsArray != NULL)
    {
      for(int i = 0; i < (int)json_array_size(atomsArray); i++)
        {
          json_t *r = json_array_get(atomsArray, i);
          if(!json_is_object(r))
            {
              fprintf(stderr, "Found a non-object item where binary reactions should be\n");
              json_decref(root);
              return 1;
            }

          json_t *xj = json_object_get(r, "x");
          if(xj == NULL || !json_is_number(xj))
            {
              std::cout << "x must be a real" << std::endl;

              if(xj != NULL)
                json_decref(root);

              return 1;
            }

          json_t *yj = json_object_get(r, "y");
          if(yj == NULL || !json_is_number(yj))
            {
              std::cout << "y must be a real" << std::endl;

              if(yj != NULL)
                json_decref(root);

              return 1;
            }

          json_t *zj = json_object_get(r, "z");
          if(zj != NULL && !json_is_number(zj))
            {
              std::cout << "z must be a real" << std::endl;

              if(zj != NULL)
                json_decref(root);

              return 1;
            }

          json_t *typej = json_object_get(r, "type");
          if(typej == NULL || !json_is_integer(typej))
            {
              std::cout << "Type must be an integer" << std::endl;

              if(typej != NULL)
                json_decref(root);

              return 1;
            }

          double x = json_number_value(xj);
          double y = json_number_value(yj);
          double z = json_number_value(zj);
      
          int type = json_integer_value(typej);

          parts.insertParticle(x, y, z, type);
        }
    }

  // Initialize particles
  std::cout << "Initializing particles" << std::endl;

  if(atomsDistributionArray != NULL)
    {
      for(int i = 0; i < (int)json_array_size(atomsDistributionArray); i++)
        {
          json_t *r = json_array_get(atomsDistributionArray, i);

          json_t *typej = json_object_get(r, "type");
          if(typej == NULL || !json_is_integer(typej))
            {
              std::cout << "Type must be an integer asdf" << std::endl;
          
              if(typej != NULL)
                json_decref(root);
          
              return 1;
            }

          json_t *numberj = json_object_get(r, "number");
          if(numberj == NULL || !json_is_integer(numberj))
            {
              std::cout << "Number must be an integer" << std::endl;
          
              if(numberj != NULL)
                json_decref(root);
          
              return 1;
            }

          int type = json_integer_value(typej);
          int number = json_integer_value(numberj);
          
          std::cout << "inserting " << number << " particles of type " << type << std::endl;
          for(int n = 0; n < number; n++)
            {
              int tt;
              for(tt = 0; tt < maxTries; tt++)
                {
                  double x = X * uniform(generator),
                    y = Y * uniform(generator),
                    z = (twoD) ? Z * 0.5 : Z * uniform(generator);

                  Particles::indexListT indexList = parts.collide(x, y, z, type);

                  if(indexList.size() < 1)
                    {
                      parts.insertParticle(x, y, z, type);

                      break;
                    }
                }

              if(tt == maxTries)
                {
                  std::cout << "Failed to insert required number of particles" << std::endl;

                  return 1;
                }
            }
        }
    }

  std::cout << "Total number of particles at time 0: " << parts.particles.size() << std::endl;

///////////////////////////////////////////
  int ns = (int)json_array_size(aTypes);
  int nspecies[ns];
//////////////////////////////////////////

  json_decref(root);

  std::ofstream output;
  std::ofstream output2; 

  output.open(outFileName);
  output2.open("pop-"+outFileName);
  //output << "reaction fired" << std::endl;

  std::cout << "Reactive Brownian Dynamics starts ... (Random Seed = " << argv[3] << ")" << std::endl;

  for(int s = 0; s < steps; s++)
    {

     ///GHETTOCODE///////////////////////////////////////////////////////
      for (int i=0;i<ns;i++) nspecies[i] = 0;
        for(auto it = parts.particles.begin(); it != parts.particles.end(); it++)
          {
              Particle &p = it->second;
              nspecies[(int)(p.type-1)] += 1;
	  }
      
      for (int i=0;i<ns;i++) output2 << nspecies[i] << " ";
      output2 << std::endl;

      if (!(s%100)) {
        std::cout << s  << " " << dt*s  << " ";
        for (int i=0;i<ns;i++) std::cout << nspecies[i] << " ";
        std::cout << std::endl;
      }
     ///GHETTOCODE///////////////////////////////////////////////////////



      //if(parts.particles.size() < 2)
      //  output << "wtf has happened?\n" << std::endl;
      if(s % printSteps == 0)
        {
          if(outputType.compare("pizza") == 0)
            {
              output << "ITEM: TIMESTEP" << std::endl;
              output << s << std::endl;

              output << "ITEM: NUMBER OF ATOMS" << std::endl;
              output << parts.particles.size() << std::endl;

              output << "ITEM: BOX BOUNDS" << std::endl;
              output << 0.0 << " " << X << std::endl;
              output << 0.0 << " " << Y << std::endl;
              output << 0.0 << " " << Z << std::endl;
              
              output << "ITEM: ATOMS" << std::endl;
            }
          else
            {
              output << parts.particles.size() << std::endl;
              output << "go spurs" << std::endl;
            }

          for(auto it = parts.particles.begin(); it != parts.particles.end(); it++)
            {
              Particle &p = it->second;
              
              if(outputType.compare("pizza") == 0)
                output << std::setprecision(17) << it->first << " " << p.type << " " << p.x << " " << p.y << " " << p.z << std::endl;
              else
                output << std::setprecision(17) << p.type << " " << p.x << " " << p.y << " " << p.z << std::endl;
            }
        }

      std::vector<int> oldParticles;
      std::set<int> deleted;
      
      for(auto it = parts.particles.begin(); it != parts.particles.end(); it++)
        {
          oldParticles.push_back(it->first);
        }

      std::random_shuffle(oldParticles.begin(), oldParticles.end());

      //Handle creation reactions
      std::pair<Reactions::mrmT::iterator, Reactions::mrmT::iterator> se = reacs.mrm.equal_range(-1);

      Reactions::mrmT::iterator it = select<Reactions::mrmT>(se.first, se.second, reacs.mrm.end(), dt, 1.0, uniform, generator);

      // If reaction found, try to make it happen!
      if(it != reacs.mrm.end())
        {
          MonatomicReaction reaction = it->second;

          //std::cout << reaction.str() << std::endl;

          int Btype = reaction.B,
            Ctype = reaction.C;

          //If there is no Btype, better be no Ctype. This is a destruction reac!
          if(Btype > 0)
            {
              int tries = 0;

              double xx, yy, zz, nx, ny, nz;
              bool initialized = false; // This is just for sanity's sake. It's possible we used xx, yy, and zz unitialized. Let's make sure this doesn't happend

              for(tries = 0; tries < maxTries; tries++)
                {
                  nx = X * uniform(generator);
                  ny = Y * uniform(generator);
                  nz = (twoD) ? Z * 0.5 : Z * uniform(generator);

                  //Check if it will be possible to insert a B atom
                  auto touching = parts.collide(nx, ny, nz, Btype);

                  // The newly inserted atom might hit the two atoms we're about to delete
                  //   But if it hits a third, or it hits an atom other than the one we're possibly
                  //     removing, reject
                  //std::cout << touching.size() << std::endl;
                  if(touching.size() >= 1)
                    {
                      if(debug)
                        printf("Failed to insert new atom, atom still moved \n");
                      
                      continue;
                    }
          
                  //Check if it will be possible to insert a C atom
                  if(Ctype > 0)
                    {
                      double r = reacs.pseps[mpp(Btype, Ctype)].sample(uniform(generator));

                      //std::cout << r << std::endl;
                  
                      double theta = (twoD) ? pi / 2.0 : pi * uniform(generator);
                      double phi = 2 * pi * uniform(generator);
                  
                      initialized = true; // We initialized xx, yy, zz
                      xx = nx + r * sin(theta) * cos(phi);
                      yy = ny + r * sin(theta) * sin(phi);
                      zz = nz + r * cos(theta);

                      //printf("%f %f %f -> %f %f %f, %f\n", x[atomi][0], x[atomi][1], x[atomi][2], xx, yy, zz, r);
                      auto touching = parts.collide(xx, yy, zz, Ctype);

                      if(touching.size() >= 1)
                        {
                          if(debug)
                            printf("Failed to insert new atom, atom still moved \n");
                          
                          continue;
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
                  // Insert the two products
                  if(Btype > 0)
                    {
                      parts.insertParticle(nx, ny, nz, Btype);
                    }

                  if(Ctype > 0)
                    {
                      if(initialized == false)
                        {
                          std::cout << "xx, yy, and zz are not initialized!" << std::endl;
                          exit(-1);
                        }

                      parts.insertParticle(xx, yy, zz, Ctype);
                    }
                }
            }
        }


      for(auto it = oldParticles.begin(); it != oldParticles.end(); it++)
        {
          // Don't try to work with an atom that has already been deleted
          if(deleted.find(*it) != deleted.end())
            continue;

          int pid = *it;

          //printf("working with particle %d\n", pid);

          Particle &p = parts.particles[pid];

          double dx = gaussian(generator) * sqrt(2 * reacs.bam[p.type].D * dt);
          double dy = gaussian(generator) * sqrt(2 * reacs.bam[p.type].D * dt);
          double dz = (twoD) ? 0.0 : gaussian(generator) * sqrt(2 * reacs.bam[p.type].D * dt);

          double nx = p.x + dx, ny = p.y + dy, nz = p.z + dz;

          //vector of ints pointing to other particles
          auto touching = parts.collide(nx, ny, nz, p.type, pid);

          // Reject the situation where three particles meet
          if(touching.size() >= 2)
            {
              //std::cout << "heyhey" << std::endl;
              //parts.move(pid, p.x + dx, p.y + dy, p.z + dz);
              
              continue;
            }
          else if(touching.size() == 0) // If there are no collisions, check if this particle can react in any way
            // Any way being: A -> B
            // or   A -> B + C
            {
              std::pair<Reactions::mrmT::iterator, Reactions::mrmT::iterator> se = reacs.mrm.equal_range(p.type);

              Reactions::mrmT::iterator it = select<Reactions::mrmT>(se.first, se.second, reacs.mrm.end(), dt, 1.0, uniform, generator);

              // If no reactions are found for this particle, accept the move
              if(it == reacs.mrm.end())
                {
                  parts.move(pid, nx, ny, nz);

                  continue;
                }

              MonatomicReaction reaction = it->second;

              //std::cout << reaction.str() << std::endl;

              int Btype = reaction.B,
                Ctype = reaction.C;

              double xx, yy, zz;
              bool initialized = false; // This is just for sanity's sake. It's possible we used xx, yy, and zz unitialized. Let's make sure this doesn't happend

              //If there is no Btype, better be no Ctype. This is a destruction reac!
              if(Btype > 0)
                {
                  int tries = 0;

                  for(tries = 0; tries < maxTries; tries++)
                    {
                      //Check if it will be possible to insert a B atom
                      auto touching = parts.collide(nx, ny, nz, Btype);

                      // The newly inserted atom might hit the two atoms we're about to delete
                      //   But if it hits a third, or it hits an atom other than the one we're possibly
                      //     removing, reject
                      //std::cout << touching.size() << std::endl;
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
                  
                          double theta = (twoD) ? pi * 0.5 : pi * uniform(generator);
                          double phi = 2 * pi * uniform(generator);
                  
                          initialized = true; // We initialized xx, yy, zz
                          xx = nx + r * sin(theta) * cos(phi);
                          yy = ny + r * sin(theta) * sin(phi);
                          zz = nz + r * cos(theta);

                          //printf("%f %f %f -> %f %f %f, %f\n", x[atomi][0], x[atomi][1], x[atomi][2], xx, yy, zz, r);
                          auto touching = parts.collide(xx, yy, zz, Ctype);

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
                }
              
              // Get rid of the reacted particle
              deleted.insert(pid);

              //std::cout << "Monatomic reaction" << pid << std::endl;
              parts.deleteParticle(pid);
                  
              // Insert the two products
              if(Btype > 0)
                {
                  parts.insertParticle(nx, ny, nz, Btype);
                }
              
              if(Ctype > 0)
                {
                  if(initialized == false)
                    {
                      std::cout << "xx, yy, and zz are not initialized!" << std::endl;
                      exit(-1);
                    }
                  
                  parts.insertParticle(xx, yy, zz, Ctype);
                }
            }
          else // Handle binary reactions
            {
              // If we're touching the wall, reject
              if(touching[0] == -1)
                {
                  continue;
                }

              // Remember that findTouching call way above? Let's get the atom id out of that
              int pjd = touching[0];

              //printf("collision\n");
              Particle &p2 = parts.particles[pjd];

              // We have an i touching a j, so build the appropriate particle pair
              ParticlePair pp = mpp(p.type, parts.particles[pjd].type);

              //std::cout << "Reacting " << p.type << " " << parts.particles[pjd].type << std::endl;
              // We could either have a regular binary reaction or an enzymatic reaction
              std::pair<Reactions::brmT::iterator, Reactions::brmT::iterator> se = reacs.brm.equal_range(pp);

              //std::cout << std::distance(se.first, se.second) << " options" << std::endl;
              //std::cout << "getout" << std::endl;

              Reactions::brmT::iterator it = select<Reactions::brmT>(se.first, se.second, reacs.brm.end(), dt, reacs.paccs[pp], uniform, generator);

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
              int Dtype = reaction.D;

              int tries = 0;
              double xx, yy, zz;
              bool initialized = false; // This is just for sanity's sake. It's possible we used xx, yy, and zz unitialized. Let's make sure this doesn't happend
              //If there is no Ctype, better be no Dtype. This is a destruction reac!

              if(Ctype > 0)
                {
                  for(tries = 0; tries < maxTries; tries++)
                    {
                      //Check if it will be possible to insert a C atom
                      auto touching = parts.collide(nx, ny, nz, Ctype);

                      // The newly inserted atom might hit the two atoms we're about to delete
                      //   But if it hits a third, or it hits an atom other than the one we're possibly
                      //     removing, reject
                      //std::cout << touching.size() << std::endl;
          
                      std::set<int> tmp = { pid, pjd };
                      if(touching.size() > 2 ||
                         !std::includes(tmp.begin(), tmp.end(), touching.begin(), touching.end()))
                        {
                          if(debug)
                            printf("Failed to insert new atom, atom still moved \n");
              
                          continue;
                        }
          
                      //Check if it will be possible to insert a D atom
                      if(Dtype > 0)
                        {
                          double r = reacs.pseps[mpp(Ctype, Dtype)].sample(uniform(generator));

                          //std::cout << r << std::endl;
                  
                          double theta = (twoD) ? pi * 0.5 : pi * uniform(generator);
                          double phi = 2 * pi * uniform(generator);
                  
                          initialized = true; // We initialized xx, yy, zz
                          xx = nx + r * sin(theta) * cos(phi);
                          yy = ny + r * sin(theta) * sin(phi);
                          zz = nz + r * cos(theta);

                          //printf("%f %f %f -> %f %f %f, %f\n", x[atomi][0], x[atomi][1], x[atomi][2], xx, yy, zz, r);
                          auto touching = parts.collide(xx, yy, zz, Dtype);

                          std::set<int> tmp = { pid, pjd };
                          if(touching.size() > 2 ||
                             !std::includes(tmp.begin(), tmp.end(), touching.begin(), touching.end()))
                            {
                              if(debug)
                                printf("Failed to insert new atom, atom still moved \n");
                              
                              continue;
                            }
                        }

                      // If we got this far, any necessary reaction products have been successfully inserted
                      break;
                    }
                }

              // If we failed to insert the new atoms, reject the move and go on
              if(tries == maxTries)
                {
                  //std::cout << "things aren't fitting" << std::endl;
                  
                  continue;
                }
              else
                { 
                  // Get rid of the reacted particle
                  deleted.insert(pid);
                  deleted.insert(pjd);
                  
                  //std::cout << "Monatomic reaction" << pid << std::endl;
                  parts.deleteParticle(pid);
                  parts.deleteParticle(pjd);
                  
                  // Insert the two products
                  if(Ctype > 0)
                    {
                      parts.insertParticle(nx, ny, nz, Ctype);
                      
                      if(Dtype > 0)
                        {
                          if(initialized == false)
                            {
                              std::cout << "xx, yy, and zz are not initialized!" << std::endl;
                              exit(-1);
                            }
                          else
                            {
                              parts.insertParticle(xx, yy, zz, Dtype);
                            }
                        }
                    }
                }

              /*if(Ctype > 0)
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
                    }*/
          
              // Time to delete the reactants
              //if(debug)
              //  printf("Trying to delete\n");
            }
        }
    }
  output.close();
  output2.close();
}
