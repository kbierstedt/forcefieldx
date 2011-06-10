/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.potential;

import static java.lang.Math.exp;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.sqrt;
import static java.lang.String.format;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Logger;

import edu.rit.pj.ParallelTeam;
import ffx.crystal.Crystal;
import ffx.crystal.ReplicatesCrystal;
import ffx.numerics.Potential;
import static ffx.numerics.VectorMath.*;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.OutOfPlaneBend;
import ffx.potential.bonded.PiOrbitalTorsion;
import ffx.potential.bonded.ROLS;
import ffx.potential.bonded.StretchBend;
import ffx.potential.bonded.Torsion;
import ffx.potential.bonded.TorsionTorsion;
import ffx.potential.bonded.UreyBradley;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldInteger;
import ffx.potential.parameters.ForceField.ForceFieldString;

/**
 * Compute the potential energy and derivatives of an AMOEBA system.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class ForceFieldEnergy implements Potential, LambdaInterface {

    private static final Logger logger = Logger.getLogger(ForceFieldEnergy.class.getName());
    private final Atom[] atoms;
    private final Crystal crystal;
    private final ParallelTeam parallelTeam;
    private final Bond bonds[];
    private final Angle angles[];
    private final StretchBend stretchBends[];
    private final UreyBradley ureyBradleys[];
    private final OutOfPlaneBend outOfPlaneBends[];
    private final Torsion torsions[];
    private final PiOrbitalTorsion piOrbitalTorsions[];
    private final TorsionTorsion torsionTorsions[];
    private final VanDerWaals vanderWaals;
    private final ParticleMeshEwald particleMeshEwald;
    protected final int nAtoms;
    protected final int nBonds;
    protected final int nAngles;
    protected final int nStretchBends;
    protected final int nUreyBradleys;
    protected final int nOutOfPlaneBends;
    protected final int nTorsions;
    protected final int nPiOrbitalTorsions;
    protected final int nTorsionTorsions;
    protected int nVanDerWaals, nPME, nGK;
    protected final boolean bondTerm;
    protected final boolean angleTerm;
    protected final boolean stretchBendTerm;
    protected final boolean ureyBradleyTerm;
    protected final boolean outOfPlaneBendTerm;
    protected final boolean torsionTerm;
    protected final boolean piOrbitalTorsionTerm;
    protected final boolean torsionTorsionTerm;
    protected final boolean vanderWaalsTerm;
    protected final boolean multipoleTerm;
    protected final boolean polarizationTerm;
    protected final boolean generalizedKirkwoodTerm;
    protected final boolean lambdaTerm;
    protected double bondEnergy, bondRMSD;
    protected double angleEnergy, angleRMSD;
    protected double stretchBendEnergy;
    protected double ureyBradleyEnergy;
    protected double outOfPlaneBendEnergy;
    protected double torsionEnergy;
    protected double piOrbitalTorsionEnergy;
    protected double torsionTorsionEnergy;
    protected double totalBondedEnergy;
    protected double vanDerWaalsEnergy;
    protected double permanentMultipoleEnergy;
    protected double polarizationEnergy;
    protected double totalElectrostaticEnergy;
    protected double totalNonBondedEnergy;
    protected double solvationEnergy;
    protected double biasEnergy;
    protected double totalEnergy;
    protected long bondTime, angleTime, stretchBendTime, ureyBradleyTime;
    protected long outOfPlaneBendTime, torsionTime, piOrbitalTorsionTime;
    protected long torsionTorsionTime, vanDerWaalsTime, electrostaticTime;
    protected long recursionTime;
    protected long totalTime;
    protected double[] optimizationScaling = null;
    /**
     * Orthogonal Space Random Walk
     */
    private double lambda;
    private boolean doCounting = true;
    private int energyCount;
    /**
     * Gas constant (in Kcal/mole/Kelvin).
     */
    public static final double R = 1.9872066e-3;
    /**
     * The first Lambda bin is centered on 0.0 (-0.005 .. 0.005).
     * The final Lambda bin is centered on 1.0 ( 0.995 .. 1.005).
     *
     * With this scheme, the maximum of biasing Gaussians
     * is at the edges.
     */
    private int lambdaBins = 101;
    /**
     * It is useful to have an odd number of bins, so that there is
     * a bin from FL=-dFL/2 to dFL/2 so that as FL approaches zero its
     * contribution to thermodynamic integration goes to zero. Otherwise
     * a contribution of zero from a L bin can only result from equal
     * sampling of the ranges -dFL to 0 and 0 to dFL.
     */
    private int FLambdaBins = 401;
    private int recursionKernel[][];
    private int biasCutoff = 5;
    private double dL = 0.01;
    private double dL_2 = dL / 2.0;
    private double dFL = 2.0;
    private double dFL_2 = dFL / 2.0;
    private double minLambda = -0.005;
    private double minFLambda = -(dFL * FLambdaBins) / 2.0;
    private double maxFLambda = minFLambda + FLambdaBins * dFL;
    private double dEdLambda = 0.0;
    private double d2EdLambda2 = 0.0;
    private double dUdXdL[] = null;
    private double biasGaussianMag = 0.005;
    private double FLambda[];
    private static final double toSeconds = 0.000000001;

    public ForceFieldEnergy(MolecularAssembly molecularAssembly) {
        parallelTeam = new ParallelTeam();
        logger.info(format(" Constructing Force Field"));
        logger.info(format("\n SMP threads:                        %10d", parallelTeam.getThreadCount()));

        // Get a reference to the sorted atom array.
        atoms = molecularAssembly.getAtomArray();
        nAtoms = atoms.length;

        ForceField forceField = molecularAssembly.getForceField();
        bondTerm = forceField.getBoolean(ForceFieldBoolean.BONDTERM, true);
        angleTerm = forceField.getBoolean(ForceFieldBoolean.ANGLETERM, true);
        stretchBendTerm = forceField.getBoolean(ForceFieldBoolean.STRBNDTERM, true);
        ureyBradleyTerm = forceField.getBoolean(ForceFieldBoolean.UREYTERM, true);
        outOfPlaneBendTerm = forceField.getBoolean(ForceFieldBoolean.OPBENDTERM, true);
        torsionTerm = forceField.getBoolean(ForceFieldBoolean.TORSIONTERM, true);
        piOrbitalTorsionTerm = forceField.getBoolean(ForceFieldBoolean.PITORSTERM, true);
        torsionTorsionTerm = forceField.getBoolean(ForceFieldBoolean.TORTORTERM, true);
        vanderWaalsTerm = forceField.getBoolean(ForceFieldBoolean.VDWTERM, true);
        multipoleTerm = forceField.getBoolean(ForceFieldBoolean.MPOLETERM, true);
        polarizationTerm = forceField.getBoolean(ForceFieldBoolean.POLARIZETERM, true);
        generalizedKirkwoodTerm = forceField.getBoolean(ForceFieldBoolean.GKTERM, false);
        lambdaTerm = forceField.getBoolean(ForceFieldBoolean.LAMBDATERM, false);

        // Define the cutoff lengths.
        double vdwOff = forceField.getDouble(ForceFieldDouble.VDW_CUTOFF, 9.0);
        double ewaldOff = forceField.getDouble(ForceFieldDouble.EWALD_CUTOFF, 7.0);
        double buff = 2.0;
        double cutOff2 = 2.0 * (max(vdwOff, ewaldOff) + buff);

        // Determine the unit cell dimensions and Spacegroup
        String spacegroup;
        double a, b, c, alpha, beta, gamma;
        boolean aperiodic;
        try {
            a = forceField.getDouble(ForceFieldDouble.A_AXIS);
            aperiodic = false;
            b = forceField.getDouble(ForceFieldDouble.B_AXIS, a);
            c = forceField.getDouble(ForceFieldDouble.C_AXIS, a);
            alpha = forceField.getDouble(ForceFieldDouble.ALPHA, 90.0);
            beta = forceField.getDouble(ForceFieldDouble.BETA, 90.0);
            gamma = forceField.getDouble(ForceFieldDouble.GAMMA, 90.0);
            spacegroup = forceField.getString(ForceFieldString.SPACEGROUP, "P1");
        } catch (Exception e) {
            logger.info(" The system will be treated as aperiodic.");
            aperiodic = true;
            spacegroup = "P1";
            /**
             * Search all atom pairs to find the largest pair-wise distance.
             */
            double xr[] = new double[3];
            double maxr = 0.0;
            for (int i = 0; i < nAtoms - 1; i++) {
                double[] xi = atoms[i].getXYZ();
                for (int j = i + 1; j < nAtoms; j++) {
                    double[] xj = atoms[j].getXYZ();
                    diff(xi, xj, xr);
                    double r = r(xr);
                    if (r > maxr) {
                        maxr = r;
                    }
                }
            }
            a = 2.0 * (maxr + ewaldOff);
            b = a;
            c = a;
            alpha = 90.0;
            beta = 90.0;
            gamma = 90.0;
        }
        Crystal unitCell = new Crystal(a, b, c, alpha, beta, gamma, spacegroup);
        unitCell.setAperiodic(aperiodic);

        /**
         * Do we need a ReplicatesCrystal?
         */
        int l = 1;
        int m = 1;
        int n = 1;
        while (unitCell.a * l < cutOff2) {
            l++;
        }
        while (unitCell.b * m < cutOff2) {
            m++;
        }
        while (unitCell.c * n < cutOff2) {
            n++;
        }

        if (l * m * n > 1 && !aperiodic) {
            this.crystal = new ReplicatesCrystal(unitCell, l, m, n);
        } else {
            this.crystal = unitCell;
        }

        logger.info(crystal.toString());

        boolean rigidHydrogens = forceField.getBoolean(ForceFieldBoolean.RIGID_HYDROGENS, false);
        double rigidScale = forceField.getDouble(ForceFieldDouble.RIGID_SCALE, 10.0);

        if (rigidScale <= 1.0) {
            rigidScale = 1.0;
        }

        logger.info("\n Bonded Terms\n");
        if (rigidHydrogens && rigidScale > 1.0) {
            logger.info(format(" Rigid hydrogens:                    %10.2f", rigidScale));
        }

        // Collect, count, pack and sort bonds.
        if (bondTerm) {
            ArrayList<ROLS> bond = molecularAssembly.getBondList();
            nBonds = bond.size();
            bonds = bond.toArray(new Bond[nBonds]);
            Arrays.sort(bonds);
            logger.info(format(" Bonds:                              %10d", nBonds));
        } else {
            nBonds = 0;
            bonds = null;
        }

        // Collect, count, pack and sort angles.
        if (angleTerm) {
            ArrayList<ROLS> angle = molecularAssembly.getAngleList();
            nAngles = angle.size();
            angles = angle.toArray(new Angle[nAngles]);
            Arrays.sort(angles);
            logger.info(format(" Angles:                             %10d", nAngles));
        } else {
            nAngles = 0;
            angles = null;
        }

        // Collect, count, pack and sort stretch-bends.
        if (stretchBendTerm) {
            ArrayList<ROLS> stretchBend = molecularAssembly.getStretchBendList();
            nStretchBends = stretchBend.size();
            stretchBends = stretchBend.toArray(new StretchBend[nStretchBends]);
            Arrays.sort(stretchBends);
            logger.info(format(" Stretch-Bends:                      %10d", nStretchBends));
        } else {
            nStretchBends = 0;
            stretchBends = null;
        }

        // Collect, count, pack and sort Urey-Bradleys.
        if (ureyBradleyTerm) {
            ArrayList<ROLS> ureyBradley = molecularAssembly.getUreyBradleyList();
            nUreyBradleys = ureyBradley.size();
            ureyBradleys = ureyBradley.toArray(new UreyBradley[nUreyBradleys]);
            Arrays.sort(ureyBradleys);
            logger.info(format(" Urey-Bradleys:                      %10d", nUreyBradleys));
        } else {
            nUreyBradleys = 0;
            ureyBradleys = null;
        }

        /**
         * Set a multiplier on the force constants of bonded terms containing
         * hydrogens.
         */
        if (rigidHydrogens) {
            if (bonds != null) {
                for (Bond bond : bonds) {
                    if (bond.containsHydrogen()) {
                        bond.setRigidScale(rigidScale);
                    }
                }
            }
            if (angles != null) {
                for (Angle angle : angles) {
                    if (angle.containsHydrogen()) {
                        angle.setRigidScale(rigidScale);
                    }
                }
            }
            if (stretchBends != null) {
                for (StretchBend stretchBend : stretchBends) {
                    if (stretchBend.containsHydrogen()) {
                        stretchBend.setRigidScale(rigidScale);
                    }
                }
            }
            if (ureyBradleys != null) {
                for (UreyBradley ureyBradley : ureyBradleys) {
                    if (ureyBradley.containsHydrogen()) {
                        ureyBradley.setRigidScale(rigidScale);
                    }
                }
            }
        }

        // Collect, count, pack and sort out-of-plane bends.
        if (outOfPlaneBendTerm) {
            ArrayList<ROLS> outOfPlaneBend = molecularAssembly.getOutOfPlaneBendList();
            nOutOfPlaneBends = outOfPlaneBend.size();
            outOfPlaneBends = outOfPlaneBend.toArray(new OutOfPlaneBend[nOutOfPlaneBends]);
            Arrays.sort(outOfPlaneBends);
            logger.info(format(" Out-of-Plane Bends:                 %10d", nOutOfPlaneBends));
        } else {
            nOutOfPlaneBends = 0;
            outOfPlaneBends = null;
        }

        // Collect, count, pack and sort torsions.
        if (torsionTerm) {
            ArrayList<ROLS> torsion = molecularAssembly.getTorsionList();
            nTorsions = torsion.size();
            torsions = torsion.toArray(new Torsion[nTorsions]);
            logger.info(format(" Torsions:                           %10d", nTorsions));
        } else {
            nTorsions = 0;
            torsions = null;
        }

        // Collect, count, pack and sort pi-orbital torsions.
        if (piOrbitalTorsionTerm) {
            ArrayList<ROLS> piOrbitalTorsion = molecularAssembly.getPiOrbitalTorsionList();
            nPiOrbitalTorsions = piOrbitalTorsion.size();
            piOrbitalTorsions = piOrbitalTorsion.toArray(new PiOrbitalTorsion[nPiOrbitalTorsions]);
            logger.info(format(" Pi-Orbital Torsions:                %10d", nPiOrbitalTorsions));
        } else {
            nPiOrbitalTorsions = 0;
            piOrbitalTorsions = null;
        }

        // Collect, count, pack and sort torsion-torsions.
        if (torsionTorsionTerm) {
            ArrayList<ROLS> torsionTorsion = molecularAssembly.getTorsionTorsionList();
            nTorsionTorsions = torsionTorsion.size();
            torsionTorsions = torsionTorsion.toArray(new TorsionTorsion[nTorsionTorsions]);
            logger.info(format(" Torsion-Torsions:                   %10d", nTorsionTorsions));
        } else {
            nTorsionTorsions = 0;
            torsionTorsions = null;
        }

        logger.info("\n Non-Bonded Terms");

        if (vanderWaalsTerm) {
            vanderWaals = new VanDerWaals(forceField, atoms, crystal, parallelTeam);
        } else {
            vanderWaals = null;
        }

        if (multipoleTerm) {
            particleMeshEwald = new ParticleMeshEwald(forceField, atoms, crystal, parallelTeam,
                                                      vanderWaals.getNeighborLists());
        } else {
            particleMeshEwald = null;
        }

        if (lambdaTerm) {
            biasCutoff = forceField.getInteger(ForceFieldInteger.LAMBDA_BIAS_CUTOFF, 5);
            biasGaussianMag = forceField.getDouble(ForceFieldDouble.BIAS_GAUSSIAN_MAG, 0.005);
            dL = forceField.getDouble(ForceFieldDouble.LAMBDA_BIN_WIDTH, 0.01);
            
            /**
             * Require modest sampling of the lambda path. 
             */
            if (dL > 0.1) {
                dL = 0.01;
            }
            
            /**
             * Many lambda bin widths do not evenly divide into 1.0; here we
             * correct for this by computing an integer number of bins, then
             * re-setting the lambda variable appropriately. Note that we also
             * choose to have an odd number of lambda bins, so that the centers
             * of the first and last bin are at 0 and 1.
             */
            lambdaBins = (int) (1.0 / dL);
            if (lambdaBins % 2 == 0) lambdaBins++;
            dL = 1.0 / (lambdaBins - 1);
            dL_2 = dL / 2.0;
            minLambda = -dL_2;
            
            /**
             * The initial number of FLambda bins does not really matter, since
             * a larger number are automatically allocated as needed. The center
             * of the central bin is at 0.
             */
            dFL = forceField.getDouble(ForceFieldDouble.FLAMBDA_BIN_WIDTH, 2.0);
            dFL_2 = dFL / 2.0;
            FLambdaBins = 401;
            minFLambda = -(dFL * FLambdaBins) / 2.0;
            maxFLambda = minFLambda + FLambdaBins * dFL;
            
            /**
             * Allocate space for the recursion kernel that stores counts.
             */
            recursionKernel = new int[lambdaBins][FLambdaBins];
            FLambda = new double[lambdaBins];
            dUdXdL = new double[nAtoms * 3];
            energyCount = 0;
            
        }
        molecularAssembly.setPotential(this);
    }

    public double energy(boolean gradient, boolean print) {
        bondTime = 0;
        angleTime = 0;
        stretchBendTime = 0;
        ureyBradleyTime = 0;
        outOfPlaneBendTime = 0;
        torsionTime = 0;
        piOrbitalTorsionTime = 0;
        torsionTorsionTime = 0;
        vanDerWaalsTime = 0;
        electrostaticTime = 0;
        recursionTime = 0;
        totalTime = System.nanoTime();

        // Zero out the potential energy of each bonded term.
        bondEnergy = 0.0;
        angleEnergy = 0.0;
        stretchBendEnergy = 0.0;
        ureyBradleyEnergy = 0.0;
        outOfPlaneBendEnergy = 0.0;
        torsionEnergy = 0.0;
        piOrbitalTorsionEnergy = 0.0;
        torsionTorsionEnergy = 0.0;
        totalBondedEnergy = 0.0;

        // Zero out bond and angle RMSDs.
        bondRMSD = 0.0;
        angleRMSD = 0.0;

        // Zero out the potential energy of each non-bonded term.
        vanDerWaalsEnergy = 0.0;
        permanentMultipoleEnergy = 0.0;
        polarizationEnergy = 0.0;
        totalElectrostaticEnergy = 0.0;
        totalNonBondedEnergy = 0.0;

        // Zero out the solvation energy.
        solvationEnergy = 0.0;

        // Zero out the recusion kernel energy.
        biasEnergy = 0.0;

        // Zero out the total potential energy.
        totalEnergy = 0.0;

        // Zero out the Cartesian coordinate gradient for each atom.
        if (gradient) {
            for (int i = 0; i < nAtoms; i++) {
                Atom atom = atoms[i];
                atom.setXYZGradient(0.0, 0.0, 0.0);
            }
        }

        if (bondTerm) {
            bondTime = System.nanoTime();
            for (int i = 0; i < nBonds; i++) {
                Bond b = bonds[i];
                bondEnergy += b.energy(gradient);
                double value = b.getValue();
                bondRMSD += value * value;
            }
            bondRMSD = sqrt(bondRMSD / bonds.length);
            bondTime = System.nanoTime() - bondTime;
        }


        if (angleTerm) {
            angleTime = System.nanoTime();
            for (int i = 0; i < nAngles; i++) {
                Angle a = angles[i];
                angleEnergy += a.energy(gradient);
                double value = a.getValue();
                angleRMSD += value * value;
            }
            angleRMSD = sqrt(angleRMSD / angles.length);
            angleTime = System.nanoTime() - angleTime;
        }

        if (stretchBendTerm) {
            stretchBendTime = System.nanoTime();
            for (int i = 0; i < nStretchBends; i++) {
                stretchBendEnergy += stretchBends[i].energy(gradient);
            }
            stretchBendTime = System.nanoTime() - stretchBendTime;
        }

        if (ureyBradleyTerm) {
            ureyBradleyTime = System.nanoTime();
            for (int i = 0; i < nUreyBradleys; i++) {
                ureyBradleyEnergy += ureyBradleys[i].energy(gradient);
            }
            ureyBradleyTime = System.nanoTime() - ureyBradleyTime;
        }

        if (outOfPlaneBendTerm) {
            outOfPlaneBendTime = System.nanoTime();
            for (int i = 0; i < nOutOfPlaneBends; i++) {
                outOfPlaneBendEnergy += outOfPlaneBends[i].energy(gradient);
            }
            outOfPlaneBendTime = System.nanoTime() - outOfPlaneBendTime;
        }

        if (torsionTerm) {
            torsionTime = System.nanoTime();
            for (int i = 0; i < nTorsions; i++) {
                torsionEnergy += torsions[i].energy(gradient);
            }
            torsionTime = System.nanoTime() - torsionTime;
        }

        if (piOrbitalTorsionTerm) {
            piOrbitalTorsionTime = System.nanoTime();
            for (int i = 0; i < nPiOrbitalTorsions; i++) {
                piOrbitalTorsionEnergy += piOrbitalTorsions[i].energy(gradient);
            }
            piOrbitalTorsionTime = System.nanoTime() - piOrbitalTorsionTime;
            torsionTorsionTime = System.nanoTime();
        }

        if (torsionTorsionTerm) {
            for (int i = 0; i < nTorsionTorsions; i++) {
                torsionTorsionEnergy += torsionTorsions[i].energy(gradient);
            }
            torsionTorsionTime = System.nanoTime() - torsionTorsionTime;
        }

        if (vanderWaalsTerm) {
            vanDerWaalsTime = System.nanoTime();
            vanDerWaalsEnergy = vanderWaals.energy(gradient, print);
            nVanDerWaals = this.vanderWaals.getInteractions();
            vanDerWaalsTime = System.nanoTime() - vanDerWaalsTime;
        }

        if (multipoleTerm) {
            electrostaticTime = System.nanoTime();
            totalElectrostaticEnergy = particleMeshEwald.energy(gradient, print);
            permanentMultipoleEnergy = particleMeshEwald.getPermanentEnergy();
            polarizationEnergy = particleMeshEwald.getPolarizationEnergy();
            nPME = particleMeshEwald.getInteractions();

            solvationEnergy = particleMeshEwald.getGKEnergy();
            nGK = particleMeshEwald.getGKInteractions();

            electrostaticTime = System.nanoTime() - electrostaticTime;
        }

        if (lambdaTerm) {
            if (doCounting) {
                energyCount++;
            }

            dEdLambda = 0.0;
            d2EdLambda2 = 0.0;
            if (vanderWaalsTerm) {
                dEdLambda = vanderWaals.getdEdL();
                d2EdLambda2 = vanderWaals.getd2EdL2();
            }
            if (multipoleTerm) {
                dEdLambda += particleMeshEwald.getdEdL();
                d2EdLambda2 += particleMeshEwald.getd2EdL2();
            }

            checkRecursionKernelSize();

            int stateBin = (int) floor((lambda - minLambda) / dL);
            if (stateBin < 0) {
                stateBin = 0;
            }
            if (stateBin >= lambdaBins) {
                stateBin = lambdaBins - 1;
            }
            int FStateBin = (int) floor((dEdLambda - minFLambda) / dFL);
            if (FStateBin == FLambdaBins) {
                FStateBin = FLambdaBins - 1;
            }
            assert (FStateBin < FLambdaBins);
            assert (FStateBin >= 0);

            /**
             * Calculate recursion kernel G(L, dEdL) and gradient.
             */
            double dGdState = 0.0;
            double dGdFState = 0.0;
            double ls2 = 2.0 / lambdaBins * 2.0 / lambdaBins;
            double FLs2 = dFL * 2.0 * dFL * 2.0;
            for (int iL = -biasCutoff; iL <= biasCutoff; iL++) {
                int lcenter = stateBin + iL;
                double deltaL = lambda - (lcenter * dL);
                double deltaL2 = deltaL * deltaL;
                // Mirror conditions for recursion kernel counts.
                int lcount = lcenter;
                double mirrorFactor = 1.0;
                if (lcount == 0 || lcount == lambdaBins - 1) {
                    mirrorFactor = 2.0;
                } else if (lcount < 0) {
                    lcount = -lcount;
                } else if (lcount > lambdaBins - 1) {
                    // Number of bins past the last bin
                    lcount -= (lambdaBins - 1);
                    // Mirror bin
                    lcount = lambdaBins - 1 - lcount;
                }
                for (int iFL = -biasCutoff; iFL <= biasCutoff; iFL++) {
                    int FLcenter = FStateBin + iFL;
                    /**
                     * If either of the following FL edge conditions are true,
                     * then there are no counts and we continue.
                     */
                    if (FLcenter < 0 || FLcenter >= FLambdaBins) {
                        continue;
                    }
                    double deltaFL = dEdLambda - (minFLambda + FLcenter * dFL + dFL_2);
                    double deltaFL2 = deltaFL * deltaFL;
                    double weight = mirrorFactor * recursionKernel[lcount][FLcenter];
                    double e = weight * biasGaussianMag
                               * exp(-deltaL2 / (2.0 * ls2))
                               * exp(-deltaFL2 / (2.0 * FLs2));
                    biasEnergy += e;
                    dGdState -= deltaL / ls2 * e;
                    dGdFState -= deltaFL / FLs2 * e;
                }
            }

            /**
             * Lambda gradient due to recursion kernel G(L, dEdL).
             */
            dEdLambda += dGdState + dGdFState * d2EdLambda2;

            /**
             * Atomic gradient due to recursion kernel G(L, dEdL).
             */
            for (int i = 0; i < 3 * nAtoms; i++) {
                dUdXdL[i] = 0.0;
            }
            getdEdXdL(dUdXdL);
            double grad[] = new double[3];
            for (int i = 0; i < nAtoms; i++) {
                Atom atom = atoms[i];
                atom.getXYZGradient(grad);
                grad[0] += dGdFState * dUdXdL[i * 3];
                grad[1] += dGdFState * dUdXdL[i * 3 + 1];
                grad[2] += dGdFState * dUdXdL[i * 3 + 2];
                atom.setXYZGradient(grad[0], grad[1], grad[2]);
            }


            // Update free energy F(L) every ~100 steps.
            if (energyCount % 100 == 0 && doCounting) {
                updateFLambda(true);
            }

            /**
             * Compute the energy and gradient for the recursion slave at F(L)
             * using interpolation.
             */
            computeRecursionSlave();

            /**
             * Log our current state.
             */
            logger.info(String.format(" Lambda %8.6f, Bin %d, G %10.4f, dE/dLambda %10.4f",
                                      lambda, stateBin, biasEnergy, dEdLambda));
            /**
             * Meta-dynamics grid counts (every ~10 steps).
             */
            if (energyCount % 10 == 0 && doCounting) {
                recursionKernel[stateBin][FStateBin]++;
            }
        }

        totalTime = System.nanoTime() - totalTime;

        totalBondedEnergy = bondEnergy + angleEnergy + stretchBendEnergy + ureyBradleyEnergy + outOfPlaneBendEnergy + torsionEnergy + piOrbitalTorsionEnergy + torsionTorsionEnergy;
        totalNonBondedEnergy = vanDerWaalsEnergy + totalElectrostaticEnergy;
        totalEnergy = totalBondedEnergy + totalNonBondedEnergy + solvationEnergy + biasEnergy;

        if (print) {
            StringBuilder sb = new StringBuilder("\n");
            if (gradient) {
                sb.append(" Computed Potential Energy and Atomic Coordinate Gradients\n");
            } else {
                sb.append(" Computed Potential Energy\n");
            }
            sb.append(this);
            logger.info(sb.toString());
        }
        return totalEnergy;
    }

    public void doEnergyCounting(boolean doCounting) {
        this.doCounting = doCounting;
    }

    /**
     * If necessary, allocate more space.
     */
    private void checkRecursionKernelSize() {
        if (dEdLambda > maxFLambda) {
            logger.info(String.format("Current F_lambda %8.2f > maximum historgram size %8.2f.",
                                      dEdLambda, maxFLambda));

            double origDeltaG = updateFLambda(false);

            int newFStateBins = FLambdaBins;
            while (minFLambda + newFStateBins * dFL < dEdLambda) {
                newFStateBins += 100;
            }
            int newRecursionKernel[][] = new int[lambdaBins][newFStateBins];
            /**
             * We have added bins above the indeces of the current counts
             * just copy them into the new array.
             */
            for (int i = 0; i < lambdaBins; i++) {
                for (int j = 0; j < FLambdaBins; j++) {
                    newRecursionKernel[i][j] = recursionKernel[i][j];
                }
            }
            recursionKernel = newRecursionKernel;
            FLambdaBins = newFStateBins;
            maxFLambda = minFLambda + dFL * FLambdaBins;
            logger.info(String.format("New historgram %8.2f to %8.2f with %d bins.\n",
                                      minFLambda, maxFLambda, FLambdaBins));

            assert (origDeltaG == updateFLambda(false));

        }
        if (dEdLambda < minFLambda) {
            logger.info(String.format("Current F_lambda %8.2f < minimum historgram size %8.2f.",
                                      dEdLambda, minFLambda));
            int offset = 100;
            while (dEdLambda < minFLambda - offset * dFL) {
                offset += 100;
            }
            int newFStateBins = FLambdaBins + offset;
            int newRecursionKernel[][] = new int[lambdaBins][newFStateBins];
            /**
             * We have added bins below the current counts,
             * so their indeces must be increased by:
             * offset = newFLBins - FLBins
             */
            for (int i = 0; i < lambdaBins; i++) {
                for (int j = 0; j < FLambdaBins; j++) {
                    newRecursionKernel[i][j + offset] = recursionKernel[i][j];
                }
            }
            recursionKernel = newRecursionKernel;
            minFLambda = minFLambda - offset * dFL;
            FLambdaBins = newFStateBins;
            logger.info(String.format("New historgram %8.2f to %8.2f with %d bins.\n",
                                      minFLambda, maxFLambda, FLambdaBins));
        }
    }

    private double updateFLambda(boolean print) {
        double freeEnergy = 0.0;
        if (print) {
            logger.info(" Count  Lambda Bin      F_Lambda Bin   <  F_L  >     dG");
        }
        for (int iL = 0; iL < lambdaBins; iL++) {
            int ulFL = -1;
            int llFL = -1;
            // Find the smallest FL bin.
            for (int jFL = 0; jFL < FLambdaBins; jFL++) {
                int count = recursionKernel[iL][jFL];
                if (count > 0) {
                    llFL = jFL;
                    break;
                }
            }
            // Find the largest FL bin.
            for (int jFL = FLambdaBins - 1; jFL >= 0; jFL--) {
                int count = recursionKernel[iL][jFL];
                if (count > 0) {
                    ulFL = jFL;
                    break;
                }
            }

            int lambdaCount = 0;
            // The FL range that has been sampled for iL*dL to (iL+1)*dL
            double lla = minFLambda + llFL * dFL;
            double ula = minFLambda + ulFL * dFL + dFL;
            if (ulFL == -1) {
                FLambda[iL] = 0.0;
                lla = 0.0;
                ula = 0.0;
            } else {
                double sumFLambda = 0.0;
                double partitionFunction = 0.0;
                for (int jFL = llFL; jFL <= ulFL; jFL++) {
                    double a = minFLambda + jFL * dFL + dFL_2;
                    double e = exp(evaluateKernel(iL, jFL) / (R * 300.0));
                    sumFLambda += a * e;
                    partitionFunction += e;
                    lambdaCount += recursionKernel[iL][jFL];
                }
                FLambda[iL] = sumFLambda / partitionFunction;
            }

            // The first and last bins are half size.
            double delta = dL;
            if (iL == 0 || iL == lambdaBins - 1) {
                delta *= 0.5;
            }
            freeEnergy += FLambda[iL] * delta;

            if (print) {
                double llL = iL * dL - dL_2;
                double ulL = llL + dL;
                if (llL < 0.0) {
                    llL = 0.0;
                }
                if (ulL > 1.0) {
                    ulL = 1.0;
                }

                if (lambdaBins <= 100) {
                    logger.info(String.format(" %5d [%4.2f %4.2f] [%7.1f %7.1f] <%8.3f> %8.3f",
                                              lambdaCount, llL, ulL, lla, ula,
                                              FLambda[iL], freeEnergy));
                } else {
                    logger.info(String.format(" %5d [%5.3f %5.3f] [%7.1f %7.1f] <%8.3f> %8.3f",
                                              lambdaCount, llL, ulL, lla, ula,
                                              FLambda[iL], freeEnergy));
                }
            }
        }
        return freeEnergy;
    }

    private void computeRecursionSlave() {
        for (int iL0 = 0; iL0 < lambdaBins - 1; iL0++) {
            int iL1 = iL0 + 1;
            /**
             * Find bin centers and values for
             * interpolation / extrapolation points.
             */
            double L0 = iL0 * dL;
            double L1 = L0 + dL;
            double FL0 = FLambda[iL0];
            double FL1 = FLambda[iL1];
            double deltaFL = FL1 - FL0;
            /**
             * If the lambda is less than or equal to the upper limit,
             * this is the final interval. Set the upper limit to L,
             * compute the partial derivative and break.
             */
            boolean done = false;
            if (lambda <= L1) {
                done = true;
                L1 = lambda;
            }
            /**
             * Upper limit - lower limit of the integral of the
             * extrapolation / interpolation.
             */
            biasEnergy += (FL0 * L1 + deltaFL * L1 * (0.5 * L1 - L0) / dL);
            biasEnergy -= (FL0 * L0 + deltaFL * L0 * (-0.5 * L0) / dL);
            if (done) {
                /**
                 * Compute the gradient d F(L) / dL at L.
                 */
                dEdLambda += FL0 + (L1 - L0) * deltaFL / dL;
                return;
            }
        }
    }

    public double evaluateKernel(int cLambda, int cF_Lambda) {
        /**
         * Compute the value of L and FL for the
         * center of the current bin.
         */
        double vL = cLambda * dL;
        double vFL = minFLambda + cF_Lambda * dFL + dFL_2;
        /**
         * Set the variances for the Gaussian bias.
         */
        double Ls2 = 2.0 * dL * 2.0 * dL;
        double FLs2 = 2.0 * dFL * 2.0 * dFL;
        double sum = 0.0;
        for (int iL = -biasCutoff; iL <= biasCutoff; iL++) {
            int Lcenter = cLambda + iL;
            double deltaL = vL - Lcenter * dL;
            double deltaL2 = deltaL * deltaL;

            // Mirror condition for Lambda counts.
            int lcount = Lcenter;
            double mirrorFactor = 1.0;
            if (lcount == 0 || lcount == lambdaBins - 1) {
                /**
                 * The width of the first and last bins is dLambda_2,
                 * so the mirror condition is to double their counts.
                 */
                mirrorFactor = 2.0;
            } else if (lcount < 0) {
                lcount = -lcount;
            } else if (lcount > lambdaBins - 1) {
                // number of bins past the last bin.
                lcount -= (lambdaBins - 1);
                // mirror bin
                lcount = lambdaBins - 1 - lcount;
            }

            for (int jFL = -biasCutoff; jFL <= biasCutoff; jFL++) {
                int FLcenter = cF_Lambda + jFL;
                /**
                 * For FLambda outside the count matrix the weight is
                 * 0 so we continue.
                 */
                if (FLcenter < 0 || FLcenter >= FLambdaBins) {
                    continue;
                }
                double deltaFL = vFL - (minFLambda + FLcenter * dFL + dFL_2);
                double deltaFL2 = deltaFL * deltaFL;
                double weight = mirrorFactor * recursionKernel[lcount][FLcenter];
                if (weight > 0) {
                    double e = weight * biasGaussianMag * exp(-deltaL2 / (2.0 * Ls2))
                               * exp(-deltaFL2 / (2.0 * FLs2));
                    sum += e;
                }
            }
        }
        return sum;
    }

    public double getTotal() {
        return totalEnergy;
    }

    public String getPDBHeaderString() {
        energy(false, false);
        StringBuilder sb = new StringBuilder();
        sb.append("REMARK   3  CALCULATED POTENTIAL ENERGY\n");
        if (bondTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "BOND STRETCHING            : ", bondEnergy, bonds.length));
            sb.append(String.format("REMARK   3   %s %g\n",
                                    "BOND RMSD                  : ", bondRMSD));
        }
        if (angleTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "ANGLE BENDING              : ", angleEnergy, angles.length));
            sb.append(String.format("REMARK   3   %s %g\n",
                                    "ANGLE RMSD                 : ", angleRMSD));
        }
        if (stretchBendTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "STRETCH-BEND               : ", stretchBendEnergy, stretchBends.length));
        }
        if (ureyBradleyTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "UREY-BRADLEY               : ", ureyBradleyEnergy, ureyBradleys.length));
        }
        if (outOfPlaneBendTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "OUT-OF-PLANE BEND          : ", outOfPlaneBendEnergy, outOfPlaneBends.length));
        }
        if (torsionTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "TORSIONAL ANGLE            : ", torsionEnergy, torsions.length));
        }
        if (piOrbitalTorsionTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "PI-ORBITAL TORSION         : ", piOrbitalTorsionEnergy, piOrbitalTorsions.length));
        }
        if (torsionTorsionTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "TORSION-TORSION            : ", torsionTorsionEnergy, torsionTorsions.length));
        }
        if (vanderWaalsTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "VAN DER WAALS              : ", vanDerWaalsEnergy, nVanDerWaals));
        }
        if (multipoleTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "ATOMIC MULTIPOLES          : ", permanentMultipoleEnergy, nPME));
        }
        if (polarizationTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "POLARIZATION               : ", polarizationEnergy, nPME));
        }
        sb.append(String.format("REMARK   3   %s %g\n",
                                "TOTAL POTENTIAL (KCAL/MOL) : ", totalEnergy));
        int nsymm = crystal.getUnitCell().spaceGroup.getNumberOfSymOps();
        if (nsymm > 1) {
            sb.append(String.format("REMARK   3   %s %g\n",
                                    "UNIT CELL POTENTIAL        : ", totalEnergy * nsymm));
        }
        if (crystal.getUnitCell() != crystal) {
            nsymm = crystal.spaceGroup.getNumberOfSymOps();
            sb.append(String.format("REMARK   3   %s %g\n",
                                    "REPLICATES CELL POTENTIAL  : ", totalEnergy * nsymm));
        }
        sb.append("REMARK   3\n");

        return sb.toString();
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("\n");
        if (bondTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f (%8.5f)\n",
                                    "Bond Streching    ", bondEnergy, bonds.length,
                                    bondTime * toSeconds, bondRMSD));
        }
        if (angleTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f (%8.5f)\n",
                                    "Angle Bending     ", angleEnergy, angles.length,
                                    angleTime * toSeconds, angleRMSD));
        }
        if (stretchBendTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Stretch-Bend      ", stretchBendEnergy,
                                    stretchBends.length, stretchBendTime * toSeconds));
        }
        if (ureyBradleyTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Urey-Bradley      ", ureyBradleyEnergy,
                                    ureyBradleys.length, ureyBradleyTime * toSeconds));
        }
        if (outOfPlaneBendTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Out-of-Plane Bend ", outOfPlaneBendEnergy,
                                    outOfPlaneBends.length, outOfPlaneBendTime * toSeconds));
        }
        if (torsionTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Torsional Angle   ", torsionEnergy, torsions.length,
                                    torsionTime * toSeconds));
        }
        if (piOrbitalTorsionTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Pi-Orbital Torsion", piOrbitalTorsionEnergy,
                                    piOrbitalTorsions.length, piOrbitalTorsionTime * toSeconds));
        }
        if (torsionTorsionTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Torsion-Torsion   ", torsionTorsionEnergy,
                                    torsionTorsions.length, torsionTorsionTime * toSeconds));
        }
        if (vanderWaalsTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Van der Waals     ", vanDerWaalsEnergy,
                                    nVanDerWaals, vanDerWaalsTime * toSeconds));
        }
        if (multipoleTerm) {
            sb.append(String.format(" %s %16.8f %12d\n",
                                    "Atomic Multipoles ", permanentMultipoleEnergy, nPME));
        }
        if (polarizationTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Polarization      ", polarizationEnergy,
                                    nPME, electrostaticTime * toSeconds));
        }

        if (generalizedKirkwoodTerm) {
            sb.append(String.format(" %s %16.8f %12d\n",
                                    "Solvation         ", solvationEnergy, nGK));
        }

        if (lambdaTerm) {
            sb.append(String.format(" %s %16.8f\n",
                                    "Bias Energy       ", biasEnergy));
        }

        sb.append(String.format("\n %s %16.8f  %s %12.3f (sec)\n",
                                "Total Potential   ", totalEnergy, "(Kcal/mole)", totalTime * toSeconds));

        int nsymm = crystal.getUnitCell().spaceGroup.getNumberOfSymOps();
        if (nsymm > 1) {
            sb.append(String.format(" %s %16.8f\n", "Unit Cell         ",
                                    totalEnergy * nsymm));
        }
        if (crystal.getUnitCell() != crystal) {
            nsymm = crystal.spaceGroup.getNumberOfSymOps();
            sb.append(String.format(" %s %16.8f\n", "Replicates Cell   ",
                                    totalEnergy * nsymm));
        }

        return sb.toString();
    }

    public Crystal getCrystal() {
        return crystal;
    }

    @Override
    public void setLambda(double lambda) {
        if (lambda <= 1.0 && lambda >= 0.0) {
            this.lambda = lambda;
            if (vanderWaalsTerm) {
                vanderWaals.setLambda(lambda);
            }
            if (multipoleTerm) {
                particleMeshEwald.setLambda(lambda);
            }
        } else {
            String message = String.format("Lambda value %8.3f is not in the range [0..1].", lambda);
            logger.warning(message);
        }
    }

    @Override
    public void setScaling(double scaling[]) {
        if (scaling != null) {
            optimizationScaling = scaling;
        } else {
            optimizationScaling = null;
        }
    }

    @Override
    public double[] getScaling() {
        return optimizationScaling;
    }

    @Override
    public double energyAndGradient(double x[], double g[]) {
        /**
         * Unscale the coordinates.
         */
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }
        setCoordinates(x);
        double e = energy(true, false);
        getGradients(g);
        /**
         * Scale the coordinates and gradients.
         */
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] *= optimizationScaling[i];
                g[i] /= optimizationScaling[i];
            }
        }
        return e;
    }

    public void getGradients(double g[]) {
        assert (g != null);
        double grad[] = new double[3];
        int index = 0;
        for (Atom a : atoms) {
            a.getXYZGradient(grad);
            double gx = grad[0];
            double gy = grad[1];
            double gz = grad[2];
            if (gx == Double.NaN || gx == Double.NEGATIVE_INFINITY || gx == Double.POSITIVE_INFINITY
                || gy == Double.NaN || gy == Double.NEGATIVE_INFINITY || gy == Double.POSITIVE_INFINITY
                || gz == Double.NaN || gz == Double.NEGATIVE_INFINITY || gz == Double.POSITIVE_INFINITY) {
                String message = format("The gradient of atom %s is (%8.3f,%8.3f,%8.3f).",
                                        a.toString(), gx, gy, gz);
                logger.warning(message);
            }
            g[index++] = gx;
            g[index++] = gy;
            g[index++] = gz;
        }
    }

    private void setCoordinates(double coords[]) {
        assert (coords != null);
        int index = 0;
        for (Atom a : atoms) {
            double x = coords[index++];
            double y = coords[index++];
            double z = coords[index++];
            a.moveTo(x, y, z);
        }
    }

    @Override
    public double[] getCoordinates(double x[]) {
        int n = getNumberOfVariables();
        if (x == null || x.length < n) {
            x = new double[n];
        }
        double xyz[] = new double[3];
        int index = 0;
        for (Atom a : atoms) {
            a.getXYZ(xyz);
            x[index++] = xyz[0];
            x[index++] = xyz[1];
            x[index++] = xyz[2];
        }
        return x;
    }

    @Override
    public double[] getMass() {
        int n = getNumberOfVariables();
        double mass[] = new double[n];
        int i = 0;
        for (Atom a : atoms) {
            double m = a.getMass();
            mass[i++] = m;
            mass[i++] = m;
            mass[i++] = m;
        }
        return mass;
    }

    @Override
    public int getNumberOfVariables() {
        return nAtoms * 3;
    }

    @Override
    public double getdEdL() {
        return dEdLambda;
    }

    @Override
    public void getdEdXdL(double gradients[]) {
        if (multipoleTerm) {
            particleMeshEwald.getdEdXdL(gradients);
        }
        if (vanderWaalsTerm) {
            vanderWaals.getdEdXdL(gradients);
        }
    }

    @Override
    public double getLambda() {
        return lambda;
    }

    @Override
    public double getd2EdL2() {
        return d2EdLambda2;
    }
}
