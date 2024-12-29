#pragma once

#include "G4VSolid.hh"
#include "geomdefs.hh"
#include "G4SurfaceMeshCGAL.hh"
#include "G4AffineTransform.hh"

/* G4HalfSpace          FLUKA       MCNP
 * --------------------------------------
 * AARBox               RPP                  1
 * RBox                 BOX                  2
 * Sphere               SPH                  3
 * CircularCone         TRC                  4
 * Ellipsoid            ELL                  5
 * Wedge                WED/RAW              6
 * Arbitrary            ARB                  7
 *
 * Plane                PLA                  8
 * XYPlane              XYP                  9
 * XZPlane              XZP                  10
 * YZPlane              YZP                  11
 *
 * CircularCylinder     RCC                  12
 * XACircularCylinder   XCC                  13
 * YACircularCylinder   YCC                  14
 * ZACircularCylinder   ZCC                  15
 *
 * EllipticCylinder     REC                  16
 * XAEllipticalCylinder XEC                  17
 * YAEllipticalCylinder YEC                  18
 * ZAEllipticalCylinder ZEC                  19
 *
 * Quadric              QUA                  20
 */

class G4VHalfSpace {

public:
  G4VHalfSpace();
  ~G4VHalfSpace();
  virtual EInside Inside(const G4ThreeVector& p) const;
  virtual G4double Sdf(const G4ThreeVector &p) const = 0;
  G4ThreeVector Normal(const G4ThreeVector& p, G4double d) const;
  virtual std::vector<G4ThreeVector> Intersection(const G4ThreeVector& p, const G4ThreeVector &d) const = 0;

  virtual void Translate(const G4ThreeVector& t) = 0;
  virtual void Rotate(const G4RotationMatrix& r) = 0 ;
  virtual void Transform(const G4AffineTransform& a) = 0;

  virtual G4SurfaceMeshCGAL* GetSurfaceMesh() = 0;

  static void QuadraticSolve(G4double a, G4double b, G4double c, G4int &nSoln, G4double &x1, G4double &x2);
  static void CubicSolve(G4double a, G4double b, G4double c, G4double d, G4int &nSoln, G4double &x1, G4double &x2, G4double &x3);
  static void QuinticSolve(G4double a, G4double b, G4double c, G4double d, G4double e, G4int &nSoln, G4double &x1, G4double &x2, G4double &x3, G4double &x4);
  static void QuadricSolve(G4double a, G4double b, G4double c, G4double e, G4double f, G4double g, G4int &nSoln, G4double &x1, G4double &x2);
  static std::vector<G4double> PolynomialSolve(std::vector<G4double>& params);

protected:
  G4SurfaceMeshCGAL *cached_mesh = nullptr;
};