#pragma once

#include <string>

class G4LoaderFluka {
public:
  G4LoaderFluka();
  ~G4LoaderFluka();

  void Load(std::string file_name);

protected:
  // states of loader
  bool free = false;
  int itransform = -1;
  bool pp_include = true;

  // preprocessor defines

  // rototranslations

  // bodies

  // regions

  // materials?

  // material assignments

};