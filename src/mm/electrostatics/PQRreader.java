package mm.electrostatics;

import java.io.IOException;
import java.util.ArrayList;

import structure.constants.Constants;
import structure.exceptions.FileFormatException;
import structure.io.ReadFile;
import structure.math.Point3f;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.parameter.AtomType;

/**
 * Generic PQR reader class for reading in protein complex, protein, residue and
 * atom information from the ATOM and HETATM entry lines in a PQR file.
 * The only checking done concerns the allowed format and range of the numbers.
 * Chemical validity must be checked elsewhere.
 * @author Abdullah Kahraman
 * @version 0.5
 * @since 0.5
 */
public class PQRreader {
    /**
     * ReadFile object holding the content of a PQR file.
     */
    private ReadFile pqrFile;
    /**
     * path to the PQR file to be converted into Java Objects.
     */
    private String filePath = "";
    /**
     * List of AtomList object that hold all atoms of all molecule entities in a
     * PQR file.
     * @see #readAllAtoms()
     */
    private ArrayList < AtomList > allAtoms = new ArrayList < AtomList >();

    /**
     * Constructor; Reads in all ATOM and HETATM entries from a PQR file.
     * @param  fileName
     *         - Path to PQR file.
     * @throws IOException if error occurs while reading infile.
     * @throws FileFormatException if ATOM or HEATM line do not conform to
     *         PQR file format.
     */
    public PQRreader(final String fileName) throws IOException,
                                                   FileFormatException {
        this.setFileName(fileName);
        this.pqrFile = new ReadFile(fileName);
        this.readAllAtoms(this.pqrFile);
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the path to the PQR file.
     * @return String object holding the path to the PQR file.
     */
    public final String getFilePath() {
        return filePath;
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the path to the PQR file.
     * @param fileName
     *        - String object holding the path to the PQR file.
     */
    private void setFileName(final String fileName) {
        this.filePath = fileName;
    }
    //--------------------------------------------------------------------------
    /**
     * General method to read in all ATOM and HETATM lines in a PQR file.
     * @param fileContent
     *        - List of String object holding each a line of a PQR file to be
     *          read in.
     * @throws FileFormatException if ATOM or HEATM line does not conform to the
     *         PQR standards.
     * @see #parseAtom(String)
     */
    private void readAllAtoms(final ArrayList < String > fileContent)
                                                    throws FileFormatException {
           int i = 0;
           AtomList atoms = new AtomList();
           for (String line : fileContent) {
               i++;
               if (line.startsWith("ATOM  ") || line.startsWith("HETATM")) {
                   Atom atom = this.parseAtom(line);
                   atoms.add(atom);
                   // NMR files hold often more than one conformation of the
                   // same protein. For the moment it should suffice to read
                   // in only the first model.
                   if (line.startsWith("END")) {
                       if (atoms.size() > 0) {
                           this.allAtoms.add(atoms);
                           atoms = new AtomList();
                       }
                   }
               }
           }
           // if file has no atom entries than just simply return empty list of
           // atoms.
           if (atoms.size() > 0) {
               this.allAtoms.add(atoms);
           } else {
             System.err.println("WARNING: No ATOM or HETATM found in input "
                              + "file");
           }
    }
    //--------------------------------------------------------------------------
    /**
     * Reads in all information of an atom in ATOM or HETATM entry lines within
     * a PQR file. Information beyond the atomic radius column are ignored.
     * @param  line
     *         - String object holding ATOM or HETATM line text in the PQR file
     * @return AtomRadius object that holds all information of the ATOM or
     *         HETATM entry line.
     * @throws FileFormatException if ATOM or HEATM line does not conform to the
     *         PQR standards.
     */
    private Atom parseAtom(final String line) throws FileFormatException {
        Atom atom = new Atom();
        try {
            atom.setFlag(line.substring(0,6));
            atom.setSerialNumber(Integer.parseInt(line.substring(6,11).trim()));
            atom.setName(line.substring(12,16));
            atom.setAlternativeLocation(line.charAt(16));
            atom.setResidueName(line.substring(17,20));
            atom.setChainId(line.charAt(21));

            int residueNumber = Integer.parseInt(line.substring(22,26).trim());
            atom.setResidueNumber(residueNumber);

            atom.setICode(line.charAt(26));

            float x = Float.parseFloat(line.substring(30,38).trim());
            float y = Float.parseFloat(line.substring(38,46).trim());
            float z = Float.parseFloat(line.substring(46,54).trim());
            atom.setXYZ(new Point3f(x, y, z));

            atom.setPartialCharge(Float.parseFloat(
                                                 line.substring(55,62).trim()));

            atom.setTemperatureFactor(Float.parseFloat(
                                                 line.substring(63,69).trim()));

            // set type of atom
            for (AtomType atomType : AtomType.values()) {
                if (atomType.getAbbreviation().trim().equals(
                                                       atom.getName().trim())) {
                    atom.setType(atomType);
                    break;
                }
            }
        } catch (Exception e) {
            throw new FileFormatException("ERROR: " + line + "; does not seem "
                                        + "to have PQR format: "
                                        + e.getMessage()
                                        + Constants.LINE_SEPERATOR);
        }
    return atom;
    }
    //--------------------------------------------------------------------------
    /**
     * Method to return a list of all ATOM entries in a PDB file.
     * @return An AtomList object holding all atoms from a PDB file.
     */
    public final AtomList getAllAtoms() {
        AtomList selection = new AtomList();
        for (AtomList atoms : allAtoms) {
            for (Atom atom : atoms) {
                if (atom.getFlag().equals("ATOM  ")) {
                    selection.add(atom);
                }
            }
        }
        return selection;
    }
}
