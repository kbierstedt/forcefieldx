//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.algorithms.groovy.test

import com.sun.xml.bind.v2.TODO
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.optimize.Minimize
import ffx.crystal.Crystal
import ffx.crystal.CrystalPotential
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import groovy.cli.picocli.CliBuilder
import org.apache.commons.io.FilenameUtils
import org.apache.commons.math3.geometry.euclidean.threed.Rotation

import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import java.util.stream.Collectors


  /**
   * Execute the script.
   */

@Command(description = " Compare crystal packings based on intermolecular distances.", name = "ffxc test.crystalSearch")
class crystalSearch extends AlgorithmsScript {


  // Default convergence criteria.
  double eps = 1.0;

  // Default number of trials.
  int nTrials = 10;

  // Create the command line parser.
  def cli = new CliBuilder(usage: ' ffxc test.crystalSearch [options] <filename1>');

  @CommandLine.Option(names = ['-h', '--help'],
          description = 'Print this help message.')
  //cli.h(longOpt: 'help', 'Print this help message.');

  @CommandLine.Option(names = ['-e', '--eps'], paramLabel = "1.0", defaultValue = "1.0",
          description = 'RMS gradient convergence criteria')
  private double eps;
  //cli.e(longOpt: 'eps', args: 1, argName: '1.0', 'RMS gradient convergence criteria');

  @CommandLine.Option(names = ['-nt', '--nTrials'], paramLabel = "10", defaultValue = "1",
          description = 'Number of random trials')
  private int nTrials;
  //cli.t(longOpt: 'trials', args: 1, argName: '10', 'Number of random trials');

  /**
   * The final argument(s) should be two or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files",
          description = 'The atomic coordinate file in PDB or XYZ format.')
  List<String> filenames = null

  private CrystalPotential potential;
  private MolecularAssembly[] assemblies;

  // The p flag is set to Direct automatically unless user specifies otherwise.
  // this is necessary since SCF does not always converge for "bad" moves/rotations

  //@CommandLine.Option(names = ['-p', '--polarization'], paramLabel = 'Direct', defaultValue = "1",
  //        description = 'polarization model: [none / direct / mutual]' )
  //  cli.p ( longOpt : 'polarization', args : 1, argName : 'Direct',
  //  'polarization model: [none / direct / mutual]' );
  //  def options = cli.parse(args);


  // List<String> arguments = options.arguments();
//  if ( options.h || arguments == null || arguments.size ( ) < 1 ) {
//    return cli.usage();
//  }

  // Read in the first command line file.
  // String filename = arguments.get(0);



//  // Load convergence criteria.
//  if ( options.e ) {
//    eps = Double.parseDouble(options.e);
//  }
//
//  // Load number of trials
//  if ( options.t ) {
//    nTrials = Integer.parseInt(options.t);
//  }
//
//  //Set polarization scheme
//  if ( options.p ) {
//    System.setProperty("polarization", options.p);
//  }

  @Override
  crystalSearch run() {
    // PSUEDO: READ IN STRUCTURES
    if (!init()) {
      return
    }

    // Ensure file exists/create molecularAssembly molecular assembly
    if (filenames != null && filenames.size() > 0) {
      assemblies = algorithmFunctions.openAll(filenames.get(0))
      activeAssembly = assemblies[0]
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return
    } else {
      assemblies = [activeAssembly]
    }

    // Open the file.
//    systems = open(filename);
//    Crystal crystal = systems[0].getCrystal().getUnitCell();
    
    
    for (String filename: filenames) {
      logger.info("\n Searching for " + filename);
      logger.info(" RMS gradient convergence criteria: " + eps);

      MolecularAssembly molecularAssembly = assemblies[0]
      Crystal crystal = molecularAssembly.getCrystal().getUnitCell()

      ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
      Minimize minimize = new Minimize(molecularAssembly, null) // TODO: Add algorithm listener
      Atom[] atoms = molecularAssembly.getAtomArray();
      // spaceGroup = crystal.spaceGroup;
      double nSymm = crystal.spaceGroup.symOps.size()
      //spaceGroup.getNumberOfSymOps();

      //Allocate memory for Translate array
      double[] translate = new double[3];

      //Allocate memory for COM
      double[] com = new double[3];

      // Used to label trial number
      int counter = 1;

      //Avogadro #
      final double avogadro = 6.02214129e23;

      // Determine the number of atoms in molecule
      int nAtoms = atoms.length;

      //Allocate memory for mass of atoms
      //variable used in determining crystal density
      double mass = 0.0;

      // Create object "random" used to generate random doubles
      // TODO: Add better function to generate random doubles 
      Random random = new Random()

      // Allocate memory for the atom coordinates after translation, rotation, and minimization
      double[][] newAtoms = new double[nAtoms][3];

      // Allocate temporary memory for the axis coordinates for each atom
      double[] finalCoords = new double[3];

      // Allocate memory for the best atom coordinates from all conducted trials
      double[][] bestAtoms = new double[nAtoms][3];

      double[][] originalAtoms = new double[nAtoms][3];

      double[] bestCrystalParameters = new double[6];

      //Allocate memory for unit cell axis lengths and angles
      double a;
      double b;
      double c;
      double alpha;
      double beta;
      double gamma;

      //Use "a" length to set apropriate translation vectors
      double max = 0.5;
      double min = -max;

      //Allocate memory for variable used in comparision with density of crystal
      double minDensity = 0.8;
      double maxDensity = 1.3;

      //Boolean used to escape while loop
      boolean likelyDensity = false;


/**********************************************************************
 *      Translate the center of mass of the molecule to the origin
 **********************************************************************/

      // This allows for proper rotation of molecule and only needs to be conducted
      // once.

      // Assign the center of mass as an array to 0.0
      com[0] = 0.0;
      com[1] = 0.0;
      com[2] = 0.0;

      // Loop over each atom to find the average X, Y, and Z values
      for (atom in atoms) {
        com[0] = com[0] + atom.getX();
        com[1] = com[1] + atom.getY();
        com[2] = com[2] + atom.getZ();
      }

      // Divide average coordinate value by the number of atoms
      com[0] = com[0] / nAtoms;
      com[1] = com[1] / nAtoms;
      com[2] = com[2] / nAtoms;

      // Calculate the translation vector for the center of mass
      crystal.toPrimaryCell(com, translate);
      translate[0] = translate[0] - com[0];
      translate[1] = translate[1] - com[1];
      translate[2] = translate[2] - com[2];

      // Move each atom and calculate the total mass of the molecule
      for (atom in atoms) {
        atom.move(translate);
        mass += atom.getMass();
      }

      // For each trial apply a rotation and translation.
      for (int i = 0; i < nTrials; i++) {
        /****************************************************
         Randomly assign lattice parameters for crystal
         ****************************************************/
        while (!likelyDensity) {
          a = getRandomNumber(14, 8, random);
          b = getRandomNumber(14, 8, random);
          c = getRandomNumber(14, 8, random);

          //Monoclinic Space group settings require alpha and gamma = 90
          alpha = 90;
          //getRandomNumber(180, 20, random);
          beta = getRandomNumber(150, 140, random);
          gamma = 90;
          //getRandomNumber(180, 20, random);

          //logger.info("The crystal parameters are: " + a + " " + b + " " + c + " " + alpha+ " " + beta + " " + gamma + " " + nSymm)
          crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma);
          double den = density(mass, nSymm, avogadro, crystal);
          //logger.info("The crystal density is: " + den + " and mass is = " + mass);
          if (den > minDensity && den < maxDensity) {
            likelyDensity = true;
          }
        }
        logger.info("---------------------- Trial Number: " + counter + " ----------------------");
        counter++;
        logger.info("Cell parameters: " + a + ' ' + b + ' ' + c + ' ' + alpha + beta); // TODO: Fix print statement with alpha and gamma 
        // Set likelyDensity to false so that on next trial run it will generate new unit
        // cell parameters.
        likelyDensity = false;

        // Apply the proposed boundary condition.
        forceFieldEnergy.setCrystal(crystal);

        /***********************************************
         Apply rotation algorithm
         ***********************************************/

        // Generate random Quaternion. In order to be random it can not be evenly distributed along the angle due to spherical rotation pattern.
        double s = random.nextDouble();
        double sigma1 = Math.sqrt(1 - s);
        double sigma2 = Math.sqrt(s);
        double theta1 = 2 * Math.PI * random.nextDouble();
        double theta2 = 2 * Math.PI * random.nextDouble();
        double w = Math.cos(theta2) * sigma2;
        double x = Math.sin(theta1) * sigma1;
        double y = Math.cos(theta1) * sigma1;
        double z = Math.sin(theta2) * sigma2;

        // Create object for rotation of molecule based on quaternion
        Rotation rotation = new Rotation(w, x, y, z, true)

        // Apply rotation to atoms
        for (int j = 0; j < nAtoms; j++) {
          atoms[j].rotate(rotation.getMatrix());
        }

        /************************************************************************
         Generate a random TRANSLATION vector i times for the a, b, and c axis
         ************************************************************************/

        double d = getRandomNumber(max, min, random);
        double e = getRandomNumber(max, min, random);
        double f = getRandomNumber(max, min, random);

        double[] translation = [d, e, f];
        logger.info("Translation vector: " + translation);
        // Apply the translation to each atom in the molecule
        for (int k = 0; k < nAtoms; k++) {
          atoms[k].move(translation);
        }

      // Comment for testing GitHub


        // Minimize the current atom configuration
        Potential g = minimize.minimize(eps);

        // Gather the new minimized atom coordinates
        newAtoms = getAtomicCoordinates(atoms, nAtoms)

        // Determine the energy associated with the current configuration
        energy = forceFieldEnergy.energy(false, true);

        // If first trial hold onto energy regardless of value for future comparison
        if (i == 0) {
          energy2 = energy;
          for (l = 0; l < nAtoms; l++) {
            for (m = 0; m < 3; m++) {
              originalAtoms[l][m] = newAtoms[l][m];
            }
          }
        }

        // If new energy is lower than previous energy, store new energy.
        else if (energy < energy2) {
          energy2 = energy;
          for (l = 0; l < nAtoms; l++) {
            for (m = 0; m < 3; m++) {
              bestAtoms[l][m] = newAtoms[l][m];
            }
            if (bestAtoms[l] == [0, 0, 0]) {
              logger.severe("Atoms are saving incorrectly as [0,0,0]");
            }
          }
          bestCrystalParameters = [a, b, c, alpha, beta, gamma];
          //logger.info("The crystal parameter is "+bestCrystalParameters[0]);
        }
        // After running 1,000 trials coords were like 180, 93,100 so i believe the tranlations were summed.
        for (n = 0; n < nAtoms; n++) {
          finalCoords[0] = originalAtoms[n][0];
          finalCoords[1] = originalAtoms[n][1];
          finalCoords[2] = originalAtoms[n][2];
          atoms[n].moveTo(finalCoords[0], finalCoords[1], finalCoords[2]);
        }
        // If the energy is not lower than previously stored energy and its the last trial
        // then we need to retrieve the lowest energy snapshot
        if (i == nTrials - 1) {
          for (n = 0; n < nAtoms; n++) {
            finalCoords[0] = bestAtoms[n][0];
            finalCoords[1] = bestAtoms[n][1];
            finalCoords[2] = bestAtoms[n][2];
            atoms[n].moveTo(finalCoords[0], finalCoords[1], finalCoords[2]);
          }
          crystal.changeUnitCellParameters(bestCrystalParameters[0], bestCrystalParameters[1],
                  bestCrystalParameters[2], bestCrystalParameters[3], bestCrystalParameters[4],
                  bestCrystalParameters[5]);
          forceFieldEnergy.setCrystal(crystal);
          if (finalCoords == [0, 0, 0]) {
            logger.severe("Atoms are saving incorrectly as [0,0,0]");
          }
          logger.info("Final energy is " + energy2);
        }
      }

      // Save lowest energy snapshot
      String ext = FilenameUtils.getExtension(filename);
      filename = FilenameUtils.removeExtension(filename);

      if (ext.toUpperCase().contains("XYZ")) {
        saveAsXYZ(new File(filename + ".xyz"));
      } else {
        saveAsPDB(systems, new File(filename + ".pdb"));
      }
    }
    return this
  }

  @Override
  public List<Potential> getPotentials() {
    if (assemblies == null) {
      return new ArrayList<Potential>();
    }
    else {
      return Arrays.stream(assemblies).
              filter { a -> a != null }.
              map { a -> a.getPotentialEnergy() }.
              filter { e -> e != null }.
              collect(Collectors.toList());
    }
  }
  // Function to provide random doubles
  double getRandomNumber(double maxShift, double minShift, rando) {
    double rand = minShift + (maxShift - minShift) * rando.nextDouble();
    return rand;
  }

  // Function to get atomic coordinates after minimization
  double[][] getAtomicCoordinates(Atom[] x, numberOfAtoms) {
    double[] Coords = new double[3]
    double[][] atomCoords = new double[numberOfAtoms][3]
    for (int j = 0; j < numberOfAtoms; j++) {
      Coords[0] = x[j].getX();
      Coords[1] = x[j].getY();
      Coords[2] = x[j].getZ();
      for (int k = 0; k < 3; k++) {
        atomCoords[j][k] = Coords[k];
      }
    }
    return atomCoords;
  }

  // Proper documentation can be found in ffx.algoritms.barostat for the density() function.
  // May need to come out later
  double density(mass, nSymm, avogadro, unitCell) {
    double density = (mass * nSymm / avogadro) * (1.0e24 / unitCell.volume);
    // logger.info(" The variables used in density calculation are: " + mass + " " + nSymm + " " + unitCell.volume + " " + density);
    return density;
  }
}



