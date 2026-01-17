#include "G4HalfSpaceTest.hh"

G4HalfSpaceTest::G4HalfSpaceTest() :
  _file_name("g4halfspace_test.dat"),
  _type(G4HalfSpaceTest::Type::isotropic),
  _source_position(G4ThreeVector(0,0,0)),
  _ntests(1) {}

G4HalfSpaceTest::G4HalfSpaceTest(std::string &file_name,
                                 Type type,
                                 G4ThreeVector source_position,
                                 std::size_t ntests) :
  _file_name(file_name),
  _type(type),
  _source_position(source_position),
  _ntests(ntests) {}

std::size_t& G4HalfSpaceTest::current_test() {
  return _current_test;
}

G4ThreeVector G4HalfSpaceTest::position() {
  return _source_position;
}

G4ThreeVector G4HalfSpaceTest::direction() {
  if(_type == Type::xdir) {
    return G4ThreeVector(1,0,0);
  }
  else if(_type == Type::ydir) {
    return G4ThreeVector(0, 1, 0);
  }
  else if(_type == Type::zdir) {
    return G4ThreeVector(0, 0, 1);
  }

  return G4ThreeVector(0,0,0);
}

void G4HalfSpaceTest::add_data(G4HalfSpaceTestDataEntry *data) {
  _data.push_back(data);
}

void G4HalfSpaceTest::print_data() {
  for(auto d : _data) {
    G4cout << G4HalfSpaceTestDataEntry::string_from_type(d->_type) << " "
           << d->_position << " " << d->_direction << " "
           << " " << (int)d->_inside << " " << d->_distance
           << " " << d->_normal << " "
           << d->_ntrialintersections << " " << d->_nintersections << std::endl;
  }
}