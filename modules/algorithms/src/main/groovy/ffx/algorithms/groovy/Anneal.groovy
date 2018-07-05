package ffx.algorithms.groovy

import org.apache.commons.io.FilenameUtils

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.AnnealOptions
import ffx.algorithms.cli.DynamicsOptions
import ffx.potential.MolecularAssembly

import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters
import picocli.CommandLine.Command

/**
 * The Anneal script.
 * <br>
 * Usage:
 * <br>
 * ffxc Anneal [options] &lt;filename&gt;
 */
@Command(description = " Run simulated annealing on a system.", name = "ffxc Anneal")
class Anneal extends AlgorithmsScript {

    @Mixin
    DynamicsOptions dynamics

    @Mixin
    AnnealOptions anneal

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = "XYZ or PDB input files.")
    private List<String> filenames

    @Override
    Anneal run() {

        if (!init()) {
            return this
        }

        dynamics.init()

        String modelfilename
        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = algorithmFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
            modelfilename = filenames.get(0)
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
        } else {
            modelfilename = activeAssembly.getFile().getAbsolutePath();
        }

        logger.info("\n Running simulated annealing on " + modelfilename + "\n")

        ffx.algorithms.SimulatedAnnealing simulatedAnnealing = new ffx.algorithms.SimulatedAnnealing(activeAssembly,
                activeAssembly.getPotentialEnergy(), activeAssembly.getProperties(),
                algorithmListener, dynamics.thermostat, dynamics.integrator)

        simulatedAnnealing.anneal(anneal.upper, anneal.low, anneal.windows, dynamics.steps, dynamics.dt)

        String ext = FilenameUtils.getExtension(modelfilename)
        modelfilename = FilenameUtils.removeExtension(modelfilename)

        if (ext.toUpperCase().contains("XYZ")) {
            algorithmFunctions.saveAsXYZ(activeAssembly, new File(modelfilename + ".xyz"));
        } else {
            algorithmFunctions.saveAsPDB(activeAssembly, new File(modelfilename + ".pdb"));
        }

        return this
    }
}