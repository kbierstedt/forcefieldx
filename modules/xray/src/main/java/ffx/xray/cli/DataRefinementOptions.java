/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.xray.cli;

import ffx.xray.RefinementMinimize;
import picocli.CommandLine;

/**
 * Represents command line options for scripts that utilize some form of crystallographic data.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public abstract class DataRefinementOptions {

    /**
     * The refinement mode to use.
     */
    protected RefinementMinimize.RefinementMode refinementMode = RefinementMinimize.RefinementMode.COORDINATES;

    /*public void init() {
        refinementMode = RefinementMinimize.parseMode(modeString);
    }*/

    /**
     * --wA or --dataWeight The weight of the data (wA).
     */
    @CommandLine.Option(names = {"--wA", "--dataWeight"}, paramLabel = "1.0",
            description = "The weight of the real space data (wA).")
    protected double wA = 1.0;

    public double getwA() {
        return wA;
    }

    // TODO: Make DiffractionData and RealSpaceData implement the same interface.
    // TODO: Create public abstract List<ThatFileType> processData(List<String> filenames, MolecularAssembly systems[])
}