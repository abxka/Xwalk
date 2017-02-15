/*
 * (C) 2010 Abdullah Kahraman
 *
 * This software is part of the open-source project "Xwalk". You can use this
 * software under the terms of the
 * Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 * (http://creativecommons.org/licenses/by-nc-sa/3.0/).
 * This means that you
 * 1.) can copy, modify, distribute the software
 * 2.) must give credit to the author
 * 3.) must not use this work for commercial purposes
 * 4.) must license derivative works under the same or a similar license.
 */

package structure.io.pdb;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.zip.DataFormatException;

import structure.constants.Constants;
import structure.exceptions.FileFormatException;
import structure.io.GzipFileReader;
import structure.io.ReadFile;
import structure.math.Point3f;
import structure.matter.Atom;
import structure.matter.AtomList;
import structure.matter.MatterUtilities;
import structure.matter.hetgroups.SmallMolecule;
import structure.matter.parameter.AminoAcidType;
import structure.matter.parameter.AtomType;
import structure.matter.protein.AminoAcid;
import structure.matter.protein.PolyPeptide;
import structure.matter.protein.PolyPeptideList;


/**
 * Generic PDB reader class for reading in protein complex, protein, residue and
 * atom information from the ATOM and HETATM entry lines in a PDB file.
 * The only checking done concerns the allowed format and range of the numbers.
 * Chemical validity must be checked elsewhere.
 * @author Abdullah Kahraman
 * @version 0.1
 * @since 0.1
 */
public class PDBreader {
    /**
     * ReadFile object holding the content of a PDB file.
     */
    private ReadFile pdbFile;
    /**
     * path to the PDB file to be converted into Java Objects.
     */
    private String filePath = "";
    /**
     * List of AtomList object that hold all atoms of all molecule entities in a
     * PDB file.
     * @see #readAllAtoms()
     */
    private ArrayList < AtomList > allAtoms = new ArrayList < AtomList >();
    //--------------------------------------------------------------------------
    /**
     * Constructor; Reads in all ATOM and HETATM entries from a PDB file.
     * @param  fileName
     *         - Path to PDB file.
     * @throws IOException if error occurs while reading infile.
     * @throws FileFormatException if ATOM or HEATM line does not conform to the
     *         <a href="http://www.wwpdb.org/documentation/format32/sect9.html">
     *         PDB standards</a>.
     */
    public PDBreader(final String fileName) throws IOException,
                                                   FileFormatException {
        this.setFileName(fileName);
        this.pdbFile = new ReadFile(fileName);
        this.readAllAtoms(this.pdbFile);
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor. Reads in all ATOM and HETATM entries from a PDB file hold
     * within a BufferedReader object.
     * @param  bufferedReader
     *         - BufferedReader object holding the PDB file.
     * @throws IOException if an error occurs while reading the BufferedReader
     *         object.
     * @throws FileFormatException if ATOM or HEATM line does not conform to the
     *         <a href="http://www.wwpdb.org/documentation/format32/sect9.html">
     *         PDB standards</a>.
     */
    public PDBreader(final BufferedReader bufferedReader)
                                                    throws IOException,
                                                           FileFormatException {
        this.pdbFile = new ReadFile(bufferedReader);
        this.readAllAtoms(this.pdbFile);
    }
    //--------------------------------------------------------------------------
    /**
     * Constructor; Reads in all ATOM and HETATM entries from a PDB file.
     * @param  readFile
     *         - ReadFile object holding the content of a PDB file.
     * @throws FileFormatException if ATOM or HEATM line does not conform to the
     *         <a href="http://www.wwpdb.org/documentation/format32/sect9.html">
     *         PDB standards</a>.
     */
    public PDBreader(final ReadFile readFile) throws FileFormatException {
        this.pdbFile = readFile;
        this.readAllAtoms(this.pdbFile);
    }
    //--------------------------------------------------------------------------
    /**
     * Creates list of PDBreader objects from the user given input file, where
     * the input file can be a list of PDB files in a tar archive (.tar),
     * compressed by GNU zip (.gz, .tar.gz, .tgz) or simply a PDB file.
     * @param infile
     *        - String object holding the path to the input file.
     * @return List of PDBreader objects each holding the content of a single
     *         PDB file.
     * @throws IOException if input file could not be read.
     * @throws FileFormatException if ATOM or HEATM line does not conform to the
     *         <a href="http://www.wwpdb.org/documentation/format32/sect9.html">
     *         PDB standards</a>.
     * @throws DataFormatException if .gz format is unknown.
     */
    public static ArrayList < PDBreader > createPDBreaders(final String infile)
                                       throws IOException, FileFormatException,
                                              DataFormatException {

        ArrayList < PDBreader > pdbReaders = new ArrayList < PDBreader >();
        if (infile.endsWith(".tar.gz") || infile.endsWith(".tgz")) {
            GzipFileReader gzip = new GzipFileReader(infile);
            TarPDBreader tarPdb = new TarPDBreader(gzip.getGZIPInputStream());
            pdbReaders.addAll(tarPdb.getPDBreaders());

        } else if (infile.endsWith(".gz") || infile.endsWith(".gzip")) {
            GzipPDBreader gzipReader = new GzipPDBreader(infile);
            pdbReaders.add(gzipReader.getPDBreader());
        } else if (infile.endsWith(".tar")) {
            TarPDBreader tarPdb = new TarPDBreader(infile);
            pdbReaders.addAll(tarPdb.getPDBreaders());
        } else {
            pdbReaders.add(new PDBreader(infile));
        }
    return pdbReaders;
    }
    //--------------------------------------------------------------------------
    /**
     * Returns the path to the PDB file.
     * @return String object holding the path to the PDB file.
     */
    public final String getFilePath() {
        return filePath;
    }
    //--------------------------------------------------------------------------
    /**
     * Sets the path to the PDB file.
     * @param fileName
     *        - String object holding the path to the PDB file.
     */
    public final void setFileName(final String fileName) {
        this.filePath = fileName;
    }
    //--------------------------------------------------------------------------
    /**
     * General method to read in all ATOM and HETATM lines in a PDB file.
     * @param fileContent
     *        - List of String object holding each a line of a PDB file to be
     *          read in.
     * @throws FileFormatException if ATOM or HEATM line does not conform to the
     *         PDB standards at
     *         {@link http://www.wwpdb.org/documentation/format32/sect9.html}
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
     * a PDB file. Information beyond the temperature column are ignored.
     * @param  line
     *         - String object holding ATOM or HETATM line text in the PDB file
     * @return AtomRadius object that holds all information of the ATOM or
     *         HETATM entry line.
     * @throws FileFormatException if ATOM or HEATM line does not conform to the
     *         PDB standards at
     *         {@link http://www.wwpdb.org/documentation/format32/sect9.html}
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

            atom.setOccupancy(Float.parseFloat(line.substring(54,60).trim()));

            float tempFact = Float.parseFloat(line.substring(60,66).trim());
            atom.setTemperatureFactor(tempFact);

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
                                        + "to have PDB format: "
                                        + e.getMessage()
                                        + Constants.LINE_SEPERATOR);
        }
    return atom;
    }
    //--------------------------------------------------------------------------
    /**
     * Method to read in all ATOM entries in a PDB file and convert all ATOM
     * entries to AminoAcidType, Protein objects and return PDB file as
     * PolyPeptideList object.
     * @param chainIds
     *        - String object of all chain IDs to be used to build up the
     *          PolyPeptideList
     * @param alternativeLocations
     *        - String object of all alternative locations to be used to build
     *          up the PolyPeptideList
     * @return A PolyPeptideList object that consists of Protein objects, which
     *         consist themselves of AminoAcidType objects, which again consists
     *         themselves of AtomRadius objects.
     */
    public final ArrayList < PolyPeptideList > getProteinComplex(
                                               final String chainIds,
                                               final String alternativeLocations
                                                               ) {
        ArrayList < PolyPeptideList > complexes =
                                            new ArrayList < PolyPeptideList >();

        StringBuffer tmp = new StringBuffer();
        for (AminoAcidType aa : AminoAcidType.values()) {
            tmp.append("#" + aa.getThreeLetterCode() + "#");
        }
        String aaTypes = tmp.toString();

        for (AtomList atoms : this.allAtoms) {
            // First put all atoms of user requested protein chains into a
            // selection hashtable.
            Hashtable < Character, AtomList > selection =
                                       new Hashtable < Character, AtomList >();
            for (Atom atom : atoms) {
                // check if atom is part of common amino acids

                if (aaTypes.indexOf("#" + atom.getResidueName() + "#") != -1
                    &&
                    (atom.getFlag().equals("ATOM  ")
                    &&
                    (chainIds.indexOf(atom.getChainId()) != -1)
                    &&
                    (alternativeLocations.indexOf(
                                                   atom.getAlternativeLocation()
                                                 ) != -1))) {
                    if (selection.get(atom.getChainId()) == null) {
                        AtomList list = new AtomList();
                        list.add(atom);
                        selection.put(atom.getChainId(), list);
                    } else {
                        selection.get(atom.getChainId()).add(atom);
                    }
                }
            }
            PolyPeptideList complex = new PolyPeptideList();
            // Get for each protein chain all amino acids.
            int rank = 1;
            char[] chains = new char[selection.size()];
            int i = 0;
            for (Enumeration < Character > e = selection.keys();
                 e.hasMoreElements();) {
                chains[i++] = e.nextElement();
            }
            Arrays.sort(chains);

            for (char chain : chains) {
                AtomList atomList = selection.get(chain);
                ArrayList < AminoAcid > aminoAcids = this.getAllAminoAcids(
                                                                        atomList
                                                                          );
                // assign rank positions to amino acids.
                for (AminoAcid aa : aminoAcids) {
                     aa.setRank(rank++);
                }
                PolyPeptide polyPeptide = new PolyPeptide(aminoAcids);
                complex.add(polyPeptide);
            }
            complex.setName(new File(this.filePath).getName());
            complexes.add(complex);
        }
        return complexes;
    }
    //--------------------------------------------------------------------------
    /**
     * Method to read in all ATOM entries in a PDB file and convert all ATOM
     * entries to AminoAcidType, Protein objects and return PDB file as
     * PolyPeptideList object.
     * @return PolyPeptideList object that consists of Protein objects, which
     *         consist themselves of AminoAcidType objects, which again consists
     *         themselves of AtomRadius objects.
     */
    public final ArrayList < PolyPeptideList > getEntireProteinComplex() {
        // First put all atoms of user requested protein chains into a selection
        // hashtable.
        return this.getProteinComplex(
                                      Constants.ALPHANUMERIC,
                                      Constants.ALPHANUMERIC
                                     );
    }
    //--------------------------------------------------------------------------
    /**
     * Method to read in all small molecule ligands in a PDB file defined by the
     * flag HETATM.
     * @return A MolecularGroup object of single SmallMolecule objects.
     */
    public final ArrayList < ArrayList < SmallMolecule >>
                                                        getAllSmallMolecules() {
        ArrayList < ArrayList < SmallMolecule >> smallMolecules =
                                new ArrayList < ArrayList < SmallMolecule >>();
        for (AtomList atoms : this.allAtoms) {
            AtomList selection = new AtomList();

            for (Atom atom : atoms) {
                if (atom.getFlag().equals("HETATM")) {
                    selection.add(atom);
                }
            }
            smallMolecules.add(this.getAllSmallMolecules(selection));
        }

        return smallMolecules;
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
    //--------------------------------------------------------------------------
    /**
     * Extracts all amino acids as AminoAcid objects from an AtomList object.
     * An amino acid is defined as a list of atoms that have an unique
     * combination of residue number and residue name.
     * @param atomList
     *        - AtomList object holding all Atom objects of a peptide/protein
     *          structure.
     * @return An array of AminoAcid objects.
     */
    private ArrayList < AminoAcid > getAllAminoAcids(final AtomList atomList) {
        ArrayList < AminoAcid > aminoAcids = new ArrayList < AminoAcid >();
        Atom preAtom = atomList.get(0);
        AtomList residueAtoms = new AtomList();

        StringBuffer tmp = new StringBuffer();
        for (AminoAcidType aa : AminoAcidType.values()) {
            tmp.append("#" + aa.getThreeLetterCode() + "#");
        }
        String aaTypes = tmp.toString();

        for (Atom atom : atomList) {
            // check if atom is part of common amino acids
            if (!MatterUtilities.equalsResidue(atom, preAtom)) {
                AminoAcid aa = new AminoAcid(residueAtoms);
                aa.setElements();
                if (aaTypes.indexOf("#" + atom.getResidueName() + "#") != -1) {
                    aa.setRank(aminoAcids.size() + 1);
                    aminoAcids.add(aa);
                }
                preAtom = atom;
                residueAtoms = new AtomList();
                residueAtoms.add(atom);
            } else {
                residueAtoms.add(atom);
            }
        }
        // last residueAtoms should also be added as an amino acid to the
        // AminoAcid ArrayList.
        AminoAcid aa = new AminoAcid(residueAtoms);
        aa.setElements();
        if (aa.getAllAtoms().size() > 0) {
            if (aaTypes.indexOf("#" + aa.getAtom(0).getResidueName() + "#")
                !=
                -1) {
                aa.setRank(aminoAcids.size() + 1);
                aminoAcids.add(aa);
            }
        }
    return aminoAcids;
    }
    //--------------------------------------------------------------------------
    /**
     * Extracts all single molecules as Molecule objects from an AtomList
     * object. A molecule is defined as a list of atoms that have an unique
     * combination of chain Id, residue number and residue name.
     * @param  atomList
     *         - AtomList object holding all Atom objects of all small
     *           molecules.
     * @return An array of SmallMolecule objects.
     */
    private ArrayList < SmallMolecule > getAllSmallMolecules(
                                                         final AtomList atomList
                                                            ) {
        ArrayList < SmallMolecule > molecules =
                                             new ArrayList < SmallMolecule >();
        if (atomList == null || atomList.size() == 0) {
            return molecules;
        }
        Atom preAtom = atomList.get(0);
        AtomList mol = new AtomList();
        for (Atom atom : atomList) {
             if (!MatterUtilities.equalsResidue(atom, preAtom)) {
                 SmallMolecule sm = new SmallMolecule(mol);
                 sm.setElements();
                 molecules.add(sm);
                 preAtom = atom;
                 mol = new AtomList();
                 mol.add(atom);
             } else {
                 mol.add(atom);
             }
        }
        SmallMolecule sm = new SmallMolecule(mol);
        sm.setElements();
        molecules.add(sm);
    return molecules;
    }
    //--------------------------------------------------------------------------
}
