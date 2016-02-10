/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
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
package ffx.numerics;

import java.util.Arrays;
import java.util.Collection;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import static org.junit.Assert.assertEquals;

/**
 * Parameterized Test of the TensorRecursion class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@RunWith(Parameterized.class)
public class TensorRecursionTest {

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{{"Test {1.1,1.2,1.3} for order 5", 1.1e0, 1.2e0, 1.3e0, 5, 56}
        });
    }
    private final double tolerance = 1.0e-14;
    private final double r[] = new double[3];
    private final double tensors[];
    private final double noStorageTensors[];
    private final int order;
    private final int tensorCount;
    private final String info;

    public TensorRecursionTest(String info, double x, double y, double z, int order, int tensorCount) {
        this.info = info;
        this.order = order;
        r[0] = x;
        r[1] = y;
        r[2] = z;
        this.tensorCount = tensorCount;
        tensors = new double[tensorCount];
        noStorageTensors = new double[tensorCount];
    }

    /**
     * Test of tensorIndex method, of class TensorRecursion.
     */
    @Test
    public void testTensorIndex() {
        int dx = 1;
        int dy = 0;
        int dz = 0;
        int expResult = 1;
        int result = TensorRecursion.tensorIndex(dx, dy, dz, order);
        assertEquals(info, expResult, result);
    }

    /**
     * Test of tensorCount method, of class TensorRecursion.
     *
     * @since 1.0
     */
    @Test
    public void testTensorCount() {
        int result = TensorRecursion.tensorCount(order);
        assertEquals(info, tensorCount, result);
    }

    /**
     * Test of tensorRecursion and noStorageTensorRecursion methods, of class
     * TensorRecursion.
     *
     * @since 1.0
     */
    @Test
    public void testTensorRecursion() {
        TensorRecursion tensorRecursion = new TensorRecursion(order);
        tensorRecursion.tensorRecursion(r, tensors);
        tensorRecursion.noStorageTensorRecursion(r, noStorageTensors);
        for (int i = 0; i < tensorCount; i++) {
            double expect = noStorageTensors[i];
            double actual = tensors[i];
            //System.out.println(String.format("%d: %10.5f %10.5f", i, expect, actual));
            assertEquals(info + " @ " + i, expect, actual, tolerance);
        }
    }
}
