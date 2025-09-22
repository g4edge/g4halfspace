#pragma once

#include <string>
#include "G4ThreeVector.hh"
#include "geomdefs.hh"

class G4HalfSpaceTestDataEntry {
public:
  enum class Type {Inside, DistanceToIn, DistanceToInDir,
                   DistanceToOut, DistanceToOutDir, Normal};

  G4HalfSpaceTestDataEntry() {};

  G4HalfSpaceTestDataEntry(Type type,
                           G4ThreeVector position,
                           G4ThreeVector direction,
                           double distance) :
                           _type(type),
                           _position(position),
                           _direction(direction),
                           _distance(distance)
                           {}

  static std::string string_from_type(Type t) {
    static const std::unordered_map<Type, std::string> type_table {
      {Type::Inside, "Inside"},
      {Type::DistanceToIn, "DistanceToIn"},
      {Type::DistanceToInDir, "DistanceToInDir"},
      {Type::DistanceToOut, "DistanceToOut"},
      {Type::DistanceToOutDir, "DistsanceToOutDir"},
      {Type::Normal,"Normal"}
    };

    auto it = type_table.find(t);
    if (it == type_table.end())
      throw std::invalid_argument("Unknown type: " );
    return it->second;
  }

  Type _type;
  G4ThreeVector _position;
  G4ThreeVector _direction;
  G4ThreeVector _normal;
  EInside _inside;
  double _distance;
  int _ntrialintersections = -1;
  int _nintersections=-1;
};

class G4HalfSpaceTest {
public:
  enum class Type {isotropic, xdir, ydir, zdir};

  static Type type_from_string(std::string& s ) {
    static const std::unordered_map<std::string, Type> string_table {
      {"isotropic", Type::isotropic},
      {"xdir", Type::xdir},
      {"ydir",  Type::ydir},
      {"zdir",  Type::zdir}
    };

    auto it = string_table.find(s);
    if (it == string_table.end())
      throw std::invalid_argument("Unknown type: " + s);
    return it->second;
  }

  static std::string string_from_type(Type t) {
    static const std::unordered_map<Type, std::string> type_table {
      {Type::isotropic, "isotropic"},
      {Type::xdir, "xdir"},
      {Type::ydir, "ydir"},
      {Type::zdir, "zdir"}
    };

    auto it = type_table.find(t);
    if (it == type_table.end())
      throw std::invalid_argument("Unknown type: " );
    return it->second;
  }

  G4HalfSpaceTest();
  G4HalfSpaceTest(std::string &file_name,
                  Type type,
                  G4ThreeVector source_position,
                  std::size_t ntests);

  std::size_t& current_test();
  G4ThreeVector position();
  G4ThreeVector direction();

  void add_data(G4HalfSpaceTestDataEntry *data);
  void print_data();

protected:
  std::string _file_name = "";
  Type _type = Type::isotropic;
  G4ThreeVector _source_position = G4ThreeVector();
  std::size_t _ntests=0;
  std::size_t _current_test=-1;

  std::vector<G4HalfSpaceTestDataEntry*> _data;
};