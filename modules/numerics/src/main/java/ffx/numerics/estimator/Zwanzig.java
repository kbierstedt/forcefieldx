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
package ffx.numerics.estimator;

import ffx.numerics.math.FFXSummaryStatistics;
import ffx.utilities.Constants;
import org.apache.commons.math3.util.FastMath;

import java.util.Arrays;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.exp;

/**
 * The Zwanzig class implements exponential averaging/free energy
 * perturbation using the Zwanzig relationship, in either the
 * forwards or backwards direction (not both).
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class Zwanzig extends SequentialEstimator {
    private static final Logger logger = Logger.getLogger(SequentialEstimator.class.getName());
    private final int nWindows;
    private final double[] dGs;
    private final double[] uncerts;
    private final double totDG;
    private final double totUncert;
    public final Directionality directionality;

    /**
     * Estimates a free energy using the Zwanzig relationship. The temperature array can be of length 1
     * if all elements are meant to be the same temperature.
     *
     * The first dimension of the energies arrays corresponds to the lambda values/windows. The
     * second dimension (can be of uneven length) corresponds to potential energies of snapshots
     * sampled from that lambda value, calculated either at that lambda value, the lambda value below,
     * or the lambda value above. The arrays energiesLow[0] and energiesHigh[n-1] is expected to be all NaN.
     *
     * @param lambdaValues   Values of lambda dynamics was run at.
     * @param energiesLow    Potential energies of trajectory L at lambda L-dL. Ignored for forwards FEP.
     * @param energiesAt     Potential energies of trajectory L at lambda L.
     * @param energiesHigh   Potential energies of trajectory L at lambda L+dL. Ignored for backwards FEP.
     * @param temperature    Temperature each lambda window was run at (single-element indicates identical temperatures).
     * @param directionality Forwards vs. backwards FEP.
     */
    public Zwanzig(double[] lambdaValues, double[][] energiesLow, double[][] energiesAt, double[][] energiesHigh, double[] temperature, Directionality directionality) {
        super(lambdaValues, energiesLow, energiesAt, energiesHigh, temperature);
        this.directionality = directionality;
        nWindows = nTrajectories - 1;

        dGs = new double[nWindows];
        uncerts = new double[nWindows];

        boolean forwards = directionality.equals(Directionality.FORWARDS);

        double cumDG = 0;
        double cumUncert = 0;
        for (int i = 0; i < nWindows; i++) {
            int windowIndex = forwards ? 0 : 1;
            windowIndex += i;
            double[] e1 = forwards ? energiesAt[windowIndex] : energiesLow[windowIndex];
            double[] e2 = forwards ? energiesHigh[windowIndex] : energiesAt[windowIndex];
            int len = e1.length;

            if (len == 0) {
                logger.info(" Skipping " + i);
                continue;
            }

            // IMPORTANT: Use the class variable temperatures, not temperature (which may be a 1-length array).
            // ALSO IMPORTANT: The -1 factor is included in here for optimization reasons.
            double beta = -temperatures[windowIndex] * Constants.R;
            double invBeta = 1.0 / beta;

            double[] deltas = new double[len];
            double[] expDeltas = new double[len];

            for (int j = 0; j < len; j++) {
                deltas[j] = e2[j] - e1[j];
                expDeltas[j] = exp(beta * deltas[j]);
            }

            FFXSummaryStatistics deltaSummary = new FFXSummaryStatistics(deltas);
            logger.info(String.format(" Mean dE for window %d: %.5f. SD of dE: %.5f", i, deltaSummary.mean, deltaSummary.sd));
            FFXSummaryStatistics expDeltaSummary = new FFXSummaryStatistics(expDeltas);
            double dG = invBeta * FastMath.log(expDeltaSummary.mean);
            dGs[i] = dG;
            cumDG += dG;
            if (len == 1) {
                uncerts[i] = 0.0;
            } else {
                double ci = deltaSummary.confidenceInterval();
                uncerts[i] = ci;
                cumUncert += (ci * ci);
            }
        }

        totDG = cumDG;
        totUncert = Math.sqrt(cumUncert);
    }

    @Override
    public boolean isBidirectional() {
        return false;
    }

    @Override
    public double getFreeEnergy() {
        return totDG;
    }

    @Override
    public double getUncertainty() {
        return totUncert;
    }

    @Override
    public double[] getWindowEnergies() {
        return Arrays.copyOf(dGs, nWindows);
    }

    @Override
    public double[] getWindowUncertainties() {
        return Arrays.copyOf(uncerts, nWindows);
    }

    public enum Directionality {
        FORWARDS, BACKWARDS;
    }
}
