# -*- coding: utf-8;
# ------------------------------------------------------------------------------
# Copyright (C) 2020 Alexander V. Popov.
#
# This source code is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
#
# This source code is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
# ------------------------------------------------------------------------------
# DNA bending angle calculator.
# ------------------------------------------------------------------------------
from __future__ import print_function
import math
import sys
import numpy as np


# ------------------------------------------------------------------------------
# PDBFile: loads a structure in PDB format and provides easy access to its
# internal data.
# ------------------------------------------------------------------------------
class PDBFile( object ):
  # Constructor: default-initializes, or loads the PDB file (if the |filename|
  # is not empty).
  def __init__( self, filename="", skip_water=False ):
    self.residues = []
    self.current_chain = 0
    self.current_residue = 0
    self.last_sequence = 0
    if filename:
      self._LoadFromFile( filename, skip_water )

  # Loads and parses the PDB file.
  def _LoadFromFile( self, filename, skip_water=False ):
    self.current_chain = 0
    self.current_residue = -sys.maxint - 1
    self.last_sequence = -sys.maxint - 1
    file_handle = open( filename, "r" )
    lines = file_handle.readlines()
    self.residues = []
    for line in lines:
      if line[:6].rstrip() in ( "ATOM", "HETATM" ):
        self._ParseAtom( line, skip_water )
      elif line[:6].rstrip() == "TER":
        self.current_chain += 1
    file_handle.close()

  # Parses an atom description.
  def _ParseAtom( self, line, skip_water ):
    residue_name = line[17:20].strip()
    if skip_water and residue_name in ( "HOH", "WAT" ):
      return
    residue_seqn = int( line[22:26].strip() )
    atom_desc = {}
    atom_desc["name"] = line[12:16].strip()
    atom_desc["x"] = float( line[30:38].strip() )
    atom_desc["y"] = float( line[38:46].strip() )
    atom_desc["z"] = float( line[46:54].strip() )
    if residue_seqn != self.last_sequence:
      # Append an atom to the residue.
      self.last_sequence = residue_seqn
      self.current_residue = len( self.residues )
      self.residues.append( {} )
      self.residues[self.current_residue]["name"] = residue_name
      self.residues[self.current_residue]["seqn"] = residue_seqn
      self.residues[self.current_residue]["chain"] = self.current_chain
      self.residues[self.current_residue]["atoms"] = {}
    self.residues[self.current_residue]["atoms"][atom_desc["name"]] = atom_desc


# ------------------------------------------------------------------------------
# DNAStructure: creates a structure of DNA based on the PDB data provided.
# Automatically detects complementary residues.
# ------------------------------------------------------------------------------
class DNAStructure( object ):
  # Constructor: initializes the empty structure, or from the PDB file (if the
  # |pdbfile| is not none).
  def __init__( self, pdbfile=None ):
    self.residues = []
    self.current_chain = 0
    self._InitTables()
    if pdbfile:
      self._CreateFromPDB( pdbfile )

  # Builds a list of DNA curve control points.
  def GetControlPoints( self ):
    # Build a list of control points from complementary residues.
    control_points_list = []
    # Sometimes we have all the residues in a single chain, and our control
    # points may then contain duplicates. Remove them now.
    control_residue_used = [ False ] * len( self.residues )
    for residue in self.residues:
      # Ignore the second chain (we account for it when processing complementary
      # residues of the first chain. Also ignore residues without a pair.
      if residue["chain"] != 0 or residue["complementary_index"] < 0:
        continue
      if control_residue_used[residue["index"]] or \
         control_residue_used[residue["complementary_index"]]:
        continue
      # Get phosphate atoms of both residues.
      try:
        p_atom0 = residue["atoms"]["P"]
        p_atom1 = self.residues[residue["complementary_index"]]["atoms"]["P"]
      except KeyError:
        continue
      control_residue_used[residue["index"]] = True
      control_residue_used[residue["complementary_index"]] = True
      control_points_list.append( [
        ( p_atom0["x"] + p_atom1["x"] ) * 0.5,
        ( p_atom0["y"] + p_atom1["y"] ) * 0.5,
        ( p_atom0["z"] + p_atom1["z"] ) * 0.5
      ] )
    return np.array( control_points_list )

  # H-bonding description table lookup indices.
  HDT_HEAVY = 0
  HDT_HYDRO = 1
  HDT_FFCODE = 2

  # FF table lookup indices.
  FFT_RM = 0
  FFT_EM = 1
  FFT_A12 = 2
  FFT_B6 = 3

  # Initializes DNA residue H-bonding descriptions and FF parameters.
  def _InitTables( self ):
    # H-bonding description table.
    self.hbond_desc = {}
    self.hbond_desc["A"] = {}
    self.hbond_desc["A"]["pair"] = [ "T", "BRU", "5IU", "8OG" ]
    self.hbond_desc["A"]["bond"] = []
    self.hbond_desc["A"]["bond"].append( ( "N6", "H61", "" ) )
    self.hbond_desc["A"]["bond"].append( ( "N6", "H62", "" ) )
    self.hbond_desc["A"]["bond"].append( ( "N1", "", "NC" ) )
    self.hbond_desc["A"]["bond"].append( ( "N3", "", "NC" ) )
    self.hbond_desc["A"]["bond"].append( ( "N7", "", "NB" ) )
    self.hbond_desc["G"] = {}
    self.hbond_desc["G"]["pair"] = [ "C", "CX2" ]
    self.hbond_desc["G"]["bond"] = []
    self.hbond_desc["G"]["bond"].append( ( "N1", "H1", "" ) )
    self.hbond_desc["G"]["bond"].append( ( "N2", "H21", "" ) )
    self.hbond_desc["G"]["bond"].append( ( "N2", "H22", "" ) )
    self.hbond_desc["G"]["bond"].append( ( "N3", "", "NC" ) )
    self.hbond_desc["G"]["bond"].append( ( "N7", "", "NB" ) )
    self.hbond_desc["G"]["bond"].append( ( "O6", "", "O" ) )
    self.hbond_desc["T"] = {}
    self.hbond_desc["T"]["pair"] = [ "A" ]
    self.hbond_desc["T"]["bond"] = []
    self.hbond_desc["T"]["bond"].append( ( "N3", "H3", "" ) )
    self.hbond_desc["T"]["bond"].append( ( "O2", "", "O" ) )
    self.hbond_desc["T"]["bond"].append( ( "O4", "", "O" ) )
    self.hbond_desc["C"] = {}
    self.hbond_desc["C"]["pair"] = [ "G", "8OG" ]
    self.hbond_desc["C"]["bond"] = []
    self.hbond_desc["C"]["bond"].append( ( "N4", "H41", "" ) )
    self.hbond_desc["C"]["bond"].append( ( "N4", "H42", "" ) )
    self.hbond_desc["C"]["bond"].append( ( "O2", "", "O" ) )
    self.hbond_desc["C"]["bond"].append( ( "N3", "", "NC" ) )
    self.hbond_desc["BRU"] = {}
    self.hbond_desc["BRU"]["pair"] = [ "A" ]
    self.hbond_desc["BRU"]["bond"] = []
    self.hbond_desc["BRU"]["bond"].append( ( "N3", "H3", "" ) )
    self.hbond_desc["BRU"]["bond"].append( ( "O2", "", "O" ) )
    self.hbond_desc["BRU"]["bond"].append( ( "O4", "", "O" ) )
    self.hbond_desc["5IU"] = {}
    self.hbond_desc["5IU"]["pair"] = [ "A" ]
    self.hbond_desc["5IU"]["bond"] = []
    self.hbond_desc["5IU"]["bond"].append( ( "N3", "H3", "" ) )
    self.hbond_desc["5IU"]["bond"].append( ( "O2", "", "O" ) )
    self.hbond_desc["5IU"]["bond"].append( ( "O4", "", "O" ) )
    self.hbond_desc["8OG"] = {}
    self.hbond_desc["8OG"]["pair"] = [ "C", "A", "CX2" ]
    self.hbond_desc["8OG"]["bond"] = []
    self.hbond_desc["8OG"]["bond"].append( ( "N1", "H1", "" ) )
    self.hbond_desc["8OG"]["bond"].append( ( "N2", "H21", "" ) )
    self.hbond_desc["8OG"]["bond"].append( ( "N2", "H22", "" ) )
    self.hbond_desc["8OG"]["bond"].append( ( "N3", "", "NC" ) )
    self.hbond_desc["8OG"]["bond"].append( ( "N7", "", "NB" ) )
    self.hbond_desc["8OG"]["bond"].append( ( "O6", "", "O" ) )
    self.hbond_desc["8OG"]["bond"].append( ( "O8", "", "O" ) )
    self.hbond_desc["CX2"] = {}
    self.hbond_desc["CX2"]["pair"] = [ "G", "8OG" ]
    self.hbond_desc["CX2"]["bond"] = []
    self.hbond_desc["CX2"]["bond"].append( ( "N4", "H41", "" ) )
    self.hbond_desc["CX2"]["bond"].append( ( "N4", "H42", "" ) )
    self.hbond_desc["CX2"]["bond"].append( ( "O2", "", "O" ) )
    self.hbond_desc["CX2"]["bond"].append( ( "N3", "", "NC" ) )
    # FF parameters table.
    self.ff_table = {}
    self.ff_table["O"] = [ 2.05, 2.5 ]
    self.ff_table["NB"] = [ 2.15, 2.0 ]
    self.ff_table["NC"] = [ 2.15, 2.0 ]
    for ffcode in self.ff_table:
      rmin_6 = pow( self.ff_table[ffcode][self.FFT_RM], 6.0 )
      self.ff_table[ffcode].append(
        self.ff_table[ffcode][self.FFT_EM] * rmin_6 * rmin_6 )
      self.ff_table[ffcode].append(
        2.0 * self.ff_table[ffcode][self.FFT_EM] * rmin_6 )

  # Creates a DNA model from the PDB data.
  def _CreateFromPDB( self, pdbfile ):
    self.residues = []
    self.current_chain = -1
    self.pdbfile_chain = -1
    for residue in pdbfile.residues:
      if not self._IsDNAResidue( residue["name"] ):
        continue
      if residue["chain"] != self.pdbfile_chain:
        self.pdbfile_chain = residue["chain"]
        self.current_chain += 1
      dna_residue = {}
      dna_residue["index"] = len( self.residues )
      dna_residue["name"] = residue["name"]
      dna_residue["seqn"] = residue["seqn"]
      dna_residue["short"] = self._ShortNameOfResidue( residue["name"] )
      dna_residue["atoms"] = residue["atoms"]
      dna_residue["chain"] = self.current_chain
      dna_residue["complementary_index"] = -1
      self.residues.append( dna_residue )
    # Do not allow more than 2 DNA chains in the structure. We probably can't
    # (don't know how to) figure out the correct double helix.
    if self.current_chain >= 2:
      raise RuntimeError( \
          "Too many DNA chains (%i)" % ( self.current_chain + 1 ) )
    if self.residues:
      self._FindComplementaryPairs()

  # Returns whether the name of the residue matches DNA residue.
  def _ShortNameOfResidue( self, residue_name ):
    name = residue_name.strip()
    if name[0] == "D":
      name = name[1:]
    if len( name ) > 1 and name[1].isdigit():
      name = name[0]
    return name

  # Returns whether the name of the residue matches DNA residue.
  def _IsDNAResidue( self, residue_name ):
    return self._ShortNameOfResidue( residue_name ) in self.hbond_desc.keys()

  def _CalculateHBondEnergy( self, atx, ath, aty, ff_code ):
    k_hb_sigma2 = 0.018
    k_hb_cutoff = 3.5
    atx_xyz = np.array( [ atx["x"], atx["y"], atx["z"] ] )
    aty_xyz = np.array( [ aty["x"], aty["y"], aty["z"] ] )
    if ath:
      ath_xyz = np.array( [ ath["x"], ath["y"], ath["z"] ] )
    else:
      # Missing the hydrogen... Assume it will be on the line connecting X and
      # Y, at the distance of 1.0 angstrom.
      ath_xyz = aty_xyz - atx_xyz
      ath_xyz = atx_xyz + ath_xyz / np.linalg.norm( ath_xyz )
    hy_delta = aty_xyz - ath_xyz
    hy_dist = np.linalg.norm( hy_delta )
    if hy_dist > k_hb_cutoff:
      return 0.0
    hx_delta = atx_xyz - ath_xyz
    hx_dist = np.linalg.norm( hx_delta )
    angle_cos = np.dot( hy_delta, hx_delta ) / ( hy_dist * hx_dist )
    energy = math.exp( -pow( angle_cos + 1.0, 2.0 ) / k_hb_sigma2 )
    rmin = self.ff_table[ff_code][self.FFT_RM]
    emin = self.ff_table[ff_code][self.FFT_EM]
    if hy_dist <= rmin:
      energy *= -emin
    else:
      ff_a12 = self.ff_table[ff_code][self.FFT_A12]
      ff_b6 = self.ff_table[ff_code][self.FFT_B6]
      hy_dist_6 = pow( hy_dist, 6.0 )
      hy_dist_12 = hy_dist_6 * hy_dist_6
      energy *= ff_a12 / hy_dist_12 - ff_b6 / hy_dist_6
    return energy

  # This will fill out indices of complementary residues.
  def _FindComplementaryPairs( self ):
    def FindComplementaryResidue( self, desc, index ):
      best_energy = 0.0
      best_index = -1
      for other in self.residues:
        other_index = other["index"]
        # Ignore self and mismatches.
        if other_index == index or other["short"] not in desc["pair"]:
          continue
        # Ignore direct neighbours.
        if abs( other_index - index ) < 2 and \
           other["chain"] == self.residues[index]["chain"]:
          continue
        try:
          other_desc = self.hbond_desc[other["short"]]
        except KeyError:
          continue
        total_energy = 0.0
        # Process X-H..Y connections.
        for bond in desc["bond"]:
          if not bond[self.HDT_HYDRO]:
            continue
          at_x = \
            self.residues[index]["atoms"][bond[self.HDT_HEAVY]]
          try:
            at_h = \
              self.residues[index]["atoms"][bond[self.HDT_HYDRO]]
          except KeyError:
            at_h = None
          for other_bond in other_desc["bond"]:
            if other_bond[self.HDT_HYDRO]:
              continue
            at_y = \
              self.residues[other_index]["atoms"][other_bond[self.HDT_HEAVY]]
            ff_code = other_bond[self.HDT_FFCODE]
            energy = self._CalculateHBondEnergy( at_x, at_h, at_y, ff_code )
            if energy > -1.0:
              continue
            total_energy += energy
        # Process Y...H-X connections.
        for bond in desc["bond"]:
          if bond[self.HDT_HYDRO]:
            continue
          at_y = self.residues[index]["atoms"][bond[self.HDT_HEAVY]]
          ff_code = bond[self.HDT_FFCODE]
          for other_bond in other_desc["bond"]:
            if not other_bond[self.HDT_HYDRO]:
              continue
            at_x = \
              self.residues[other_index]["atoms"][other_bond[self.HDT_HEAVY]]
            try:
              at_h = \
                self.residues[other_index]["atoms"][other_bond[self.HDT_HYDRO]]
            except KeyError:
              at_h = None
            energy = self._CalculateHBondEnergy( at_x, at_h, at_y, ff_code )
            if energy > -1.0:
              continue
            total_energy += energy
        if total_energy < best_energy:
          best_energy = total_energy
          best_index = other_index
      return best_index
    for residue in self.residues:
      if residue["complementary_index"] >= 0:
        continue
      current_index = residue["index"]
      other_index = -1
      try:
        desc = self.hbond_desc[residue["short"]]
        other_index = FindComplementaryResidue( self, desc, current_index )
      except KeyError:
        pass
      if other_index >= 0:
        residue["complementary_index"] = other_index
        self.residues[other_index]["complementary_index"] = current_index

  # Print a list of complimentary residues (debug).
  def _PrintResidueList( self ):
    for residue in self.residues:
      if residue["complementary_index"] < 0:
        print( "[%i]%s-%i : no pair" % \
              ( residue["chain"], residue["name"], residue["seqn"] ) )
      else:
        print( "[%i]%s-%i : [%i]%s-%i" % \
              ( residue["chain"], residue["name"], residue["seqn"],
                self.residues[residue["complementary_index"]]["chain"],
                self.residues[residue["complementary_index"]]["name"],
                self.residues[residue["complementary_index"]]["seqn"] ) )


# ------------------------------------------------------------------------------
# PolynomialRegression: for a given number of control points, builds a curve,
# calculates derivatives and angle between tangents.
# ------------------------------------------------------------------------------
class PolynomialRegression( object ):
  # Constructor: verifies the input, assigns control points and calculates their mean.
  def __init__( self, poly_degree, control_points ):
    assert poly_degree >= 2, "Polynomial degree must be >= 2"
    assert len( control_points ) >= 6, \
        "Too few control points (%i)" % len( control_points )
    assert control_points.shape[1] == 3, "Control points must be Nx3"
    self.poly_degree = poly_degree
    self.poly_coeffs = np.zeros( poly_degree + 1 )
    self.control_points = control_points
    self.num_control_points = len( self.control_points )

  def DoRegression( self ):
    self._PutControlPointsOn2DPlane()
    self._ConstructPolynomial()

  def CalcDerivatives( self ):
    derivative0 = self._CalcDerivative( 0 )
    derivative0 += self._CalcDerivative( 1 )
    derivative0 += self._CalcDerivative( 2 )
    derivative0 /= 3.0
    derivative1 = self._CalcDerivative( self.num_control_points-3 )
    derivative1 += self._CalcDerivative( self.num_control_points-2 )
    derivative1 += self._CalcDerivative( self.num_control_points-1 )
    derivative1 /= 3.0
    return ( derivative0, derivative1 )

  def CalcAngle( self, derivatives ):
    assert len( derivatives ) == 2, "Exactly two derivatives required"
    vec1 = np.sqrt( derivatives[0] * derivatives[0] + 1.0 )
    vec2 = np.sqrt( derivatives[1] * derivatives[1] + 1.0 )
    fdot = abs( derivatives[0] * derivatives[1] + 1.0 ) / ( vec1 * vec2 )
    return np.arccos( fdot ) * 180.0 / np.pi

  def _PutControlPointsOn2DPlane( self ):
    # Project control points onto the plane using PCA.
    cov = np.cov( self.control_points, rowvar=False )
    evals, evecs = np.linalg.eig( cov )
    plane_normal = evecs[:, evals.argmin()]
    plane_normal /= np.linalg.norm( plane_normal )
    plane_dist = -np.dot( plane_normal, self.control_points.mean( axis=0 ) )
    offsets = np.dot( self.control_points, plane_normal ) + plane_dist
    self.control_points -= plane_normal * offsets.reshape( -1, 1 )
    # Now transform the projected plane, so that the z-component of control
    # points will be almost zeroed. This will also align control points in
    # such a way that the curve will run from left to right.
    origin = ( self.control_points[self.num_control_points-1] + \
               self.control_points[0] ) * 0.5
    rightdir = self.control_points[self.num_control_points-1] - \
               self.control_points[0]
    rightdir /= np.linalg.norm( rightdir )
    forwarddir = np.cross( rightdir, plane_normal )
    matrix = np.array( [ rightdir, forwarddir, plane_normal ] ).transpose()
    self.control_points = np.dot( self.control_points - origin, matrix )

  def _SumXk( self, order ):
    return np.sum( np.power( self.control_points[:, 0], order ) ) \
           if order > 0 else self.num_control_points

  def _SumXkY( self, order ):
    return np.sum( np.power( self.control_points[:, 0], order ) * \
                   self.control_points[:, 1] ) if order > 0 else \
                   np.sum( self.control_points[:, 1] )

  def _ConstructPolynomial( self ):
    numk = self.poly_degree + 1
    # Build a symmetric least-squares matrix and a B vector.
    lsm = np.zeros( ( numk, numk ) )
    vecb = np.zeros( numk )
    for i in range( 0, numk ):
      for j in range( i, numk ):
        value = lsm[i-1][j+1] if ( i > 0 and j < numk-1 ) \
                              else self._SumXk( i + j )
        lsm[i][j] = lsm[j][i] = value
      vecb[i] = self._SumXkY( i )
    # Now we have: M * polyCoeffs = B; this linear equation system must be
    # solved for the polynomial coefficients.
    self.poly_coeffs = np.linalg.solve( lsm, vecb )

  def _CalcDerivative( self, control_index ):
    # Calculate a derivative at the given control point.
    control_point = self.control_points[control_index]
    derivative = self.poly_coeffs[1]
    for k in range( 1, self.poly_degree ):
      derivative += ( k + 1 ) * self.poly_coeffs[k+1] * \
                    pow( control_point[0], k )
    return derivative

  def _DebugControlPoints( self ):
    with open( "debug.pdb", "w" ) as f_output:
      at_counter = 1
      for point in self.control_points:
        f_output.write( "%-6s%5d %-4s %-4s%c%4d    %8.3f%8.3f%8.3f\n" % \
                        ( "ATOM", at_counter, "P", "CTL", "A", 1, \
                          point[0], point[1], point[2] ) )
        at_counter += 1
      f_output.write( "%-6s%5d\n" % ( "TER", at_counter ) )


# ------------------------------------------------------------------------------
# Main function: loads a PDB file, build a DNA model, gets a list of control
# points from it, performs polynomial regression for them, and finally
# calculates the DNA bending angle. This function will print the value of the
# angle, and return zero. On error, it will return a non-zero exit code.
# ------------------------------------------------------------------------------
def Main():
  if len( sys.argv ) < 2:
    raise RuntimeError( "Input PDB file not specified." )
  # Load the PDB file.
  pdbfile = PDBFile( sys.argv[1], skip_water=True )
  # Create the DNA.
  dna = DNAStructure( pdbfile )
  # Initialize the calculator.
  calculator = PolynomialRegression( 3, dna.GetControlPoints() )
  calculator.DoRegression()
  # Calculate the angle.
  print( "%g" % calculator.CalcAngle( calculator.CalcDerivatives() ) )
  return 0


if __name__ == "__main__":
  SCRIPT_ERROR_CODE = 1
  try:
    SCRIPT_ERROR_CODE = Main()
  except ( AssertionError, IOError, RuntimeError ) as e:
    print( "ERROR: %s" % str( e ) )
  sys.exit( SCRIPT_ERROR_CODE )
