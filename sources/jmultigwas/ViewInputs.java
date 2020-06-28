package jmultigwas;

import java.awt.Component;
import java.io.File;
import java.util.prefs.Preferences;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

public class ViewInputs extends javax.swing.JPanel {

    // Attributes
    Model model;
    Controller controller;
    String outputDir;
    Preferences prefs;

    final JFileChooser fc = new JFileChooser();
    String LAST_USED_FOLDER = "";

    // Methods
    public ViewInputs(Controller controller, Model model) {
        this.controller = controller;
        this.model = model;
        initComponents();
        setDefaults();
        fieldOutputDir.setText(getLastDir().toString());
        //fieldGeno.setText(getLastDir().toString());
        //fieldPheno.setText(getLastDir().toString());

    }

    public File getLastDir() {
        prefs = Preferences.userRoot().node(getClass().getName());
        File curFile = new File(prefs.get(LAST_USED_FOLDER, new File(".").getAbsolutePath()));
        return curFile;
    }

    public void setDefaults() {
        fieldOutputDir.setText(System.getProperty("user.home"));
        fieldPloidy.setSelectedIndex(0);
        fieldGeno.setText("");
        fieldPheno.setText("");
        fieldSignificance.setText("0.05");
        fieldCorrection.setSelectedIndex(0);
        fieldModel.setSelectedIndex(0);
        fieldSNPs.setSelectedIndex(9);
        fieldFiltering.setSelectedIndex(0);
        fieldFilterMAF.setText("0.01");
        fieldFilterMIND.setText("0.1");
        fieldFilterGENO.setText("0.1");
        fieldFilterHWE.setText("1e-10");
    }

    public void clearInputs() {
        fieldPloidy.setSelectedIndex(0);
        fieldOutputDir.setText("");
        fieldGeno.setText("");
        fieldPheno.setText("");
        fieldSignificance.setText("");
        fieldCorrection.setSelectedIndex(0);
        fieldModel.setSelectedIndex(0);
        fieldSNPs.setSelectedIndex(9);
        fieldFiltering.setSelectedIndex(0);
        fieldFilterMAF.setText("");
        fieldFilterMIND.setText("");
        fieldFilterGENO.setText("");
        fieldFilterHWE.setText("");
        setEnabledInputs(false);
    }

    public void setEnabledInputs(Boolean flag) {
        Component components[] = this.getComponents();

        for (Component c : panelPaths.getComponents()) {
            c.setEnabled(flag);
        }

        for (Component c : panelParameters.getComponents()) {
            c.setEnabled(flag);
        }

        for (Component c : panelFilters.getComponents()) {
            c.setEnabled(flag);
        }

        fieldOutputDir.setEnabled(true);
        labelOutputDir.setEnabled(true);
        selOutputDir.setEnabled(true);
    }

    public String getInputValues() {
        StringBuffer txt = new StringBuffer(500);
        String ln = System.lineSeparator();
        txt.append("default:" + ln);
        txt.append("  ploidy               : " + fieldPloidy.getSelectedItem().toString() + ln);
        txt.append("  genotypeFile         : " + '"' + fieldGeno.getText() + '"' + ln);
        txt.append("  phenotypeFile        : " + '"' + fieldPheno.getText() + '"' + ln);        
        txt.append("  significanceLevel    : " + fieldSignificance.getText() + ln);
        txt.append("  correctionMethod     : " + '"' + fieldCorrection.getSelectedItem().toString() + '"' + ln);
        txt.append("  gwasModel            : " + fieldModel.getSelectedItem().toString() + ln);
        txt.append("  nBest                : " + fieldSNPs.getSelectedItem().toString() + ln);
        txt.append("  filtering            : " + fieldFiltering.getSelectedItem().toString() + ln);
        txt.append("  MAF                  : " + fieldFilterMAF.getText() + ln);
        txt.append("  MIND                 : " + fieldFilterMIND.getText() + ln);
        txt.append("  GENO                 : " + fieldFilterGENO.getText() + ln);
        txt.append("  HWE                  : " + fieldFilterHWE.getText() + ln);
        txt.append("  tools                : " + '"' + controller.getToolsToRun() + '"' + ln);

        return txt.toString();
    }

    public boolean checkCompleteInfo() {
        if (fieldOutputDir.getText().equals ("")) return false;
        if (fieldGeno.getText().equals ("")) return false;
        if (fieldPheno.getText().equals ("")) return false;
        if (fieldSignificance.getText().equals ("")) return false;
        if (fieldFilterMAF.getText().equals ("")) return false;
        if (fieldFilterMIND.getText().equals ("")) return false;
        if (fieldFilterGENO.getText().equals ("")) return false;
        if (fieldFilterHWE.getText().equals ("")) return false;
        if (controller.getToolsToRun().equals ("")) return false;
        System.out.println("Success: complete information!!");
        return true;
    }

    public void runApplication() {
        String outputDir = fieldOutputDir.getText();
        String fieldsText = getInputValues();

        //controller.onRunApplication(outputDir, fieldsText);

        //model.runApplication(outputDir, fieldsText.toString());
    }

    public String getOutputDir() {
        outputDir = fieldOutputDir.getText();
        return outputDir;
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        panelInputs = new javax.swing.JPanel();
        panelFilesTitle = new javax.swing.JPanel();
        jLabel14 = new javax.swing.JLabel();
        panelPaths = new javax.swing.JPanel();
        labelOutputDir = new javax.swing.JLabel();
        fieldOutputDir = new javax.swing.JTextField();
        fieldPheno = new javax.swing.JTextField();
        jLabel2 = new javax.swing.JLabel();
        fieldGeno = new javax.swing.JTextField();
        jLabel12 = new javax.swing.JLabel();
        selOutputDir = new javax.swing.JButton();
        selGenotypeBt = new javax.swing.JButton();
        selPhenotypeBt = new javax.swing.JButton();
        panelParametersTitle = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        panelParameters = new javax.swing.JPanel();
        jLabel8 = new javax.swing.JLabel();
        fieldModel = new javax.swing.JComboBox<>();
        jLabel3 = new javax.swing.JLabel();
        fieldSignificance = new javax.swing.JTextField();
        jLabel4 = new javax.swing.JLabel();
        fieldCorrection = new javax.swing.JComboBox<>();
        jLabel5 = new javax.swing.JLabel();
        fieldSNPs = new javax.swing.JComboBox<>();
        jLabel6 = new javax.swing.JLabel();
        fieldFiltering = new javax.swing.JComboBox<>();
        labelPloidy = new javax.swing.JLabel();
        fieldPloidy = new javax.swing.JComboBox<>();
        panelFiltersTitle = new javax.swing.JPanel();
        panelFilters = new javax.swing.JPanel();
        fieldFilterMIND = new javax.swing.JTextField();
        fieldFilterGENO = new javax.swing.JTextField();
        fieldFilterHWE = new javax.swing.JTextField();
        fieldFilterMAF = new javax.swing.JTextField();
        jLabel7 = new javax.swing.JLabel();
        jLabel15 = new javax.swing.JLabel();
        jLabel16 = new javax.swing.JLabel();
        jLabel17 = new javax.swing.JLabel();
        jLabel18 = new javax.swing.JLabel();

        setBackground(new java.awt.Color(153, 153, 255));
        setPreferredSize(new java.awt.Dimension(780, 650));
        setLayout(new java.awt.BorderLayout());

        panelInputs.setBackground(javax.swing.UIManager.getDefaults().getColor("Button.focus"));
        panelInputs.setPreferredSize(new java.awt.Dimension(780, 650));

        panelFilesTitle.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        panelFilesTitle.setLayout(new java.awt.BorderLayout());

        jLabel14.setBackground(javax.swing.UIManager.getDefaults().getColor("Button.disabledToolBarBorderBackground"));
        jLabel14.setText("Input/Output:");
        jLabel14.setOpaque(true);
        panelFilesTitle.add(jLabel14, java.awt.BorderLayout.PAGE_START);

        panelPaths.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));

        labelOutputDir.setText("Output Folder:");

        fieldOutputDir.setText("/home/lg/AAA");
        fieldOutputDir.setToolTipText("Select genotype file");

        fieldPheno.setToolTipText("Select phenotype file");
        fieldPheno.setPreferredSize(new java.awt.Dimension(90, 19));

        jLabel2.setText("Input Phenotype:");

        fieldGeno.setToolTipText("Select genotype file");
        fieldGeno.setPreferredSize(new java.awt.Dimension(90, 19));

        jLabel12.setText("Input Genotype:");

        selOutputDir.setText("Select...");
        selOutputDir.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                selOutputDirActionPerformed(evt);
            }
        });

        selGenotypeBt.setText("Select...");
        selGenotypeBt.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                selGenotypeBtActionPerformed(evt);
            }
        });

        selPhenotypeBt.setText("Select...");
        selPhenotypeBt.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                selPhenotypeBtActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout panelPathsLayout = new javax.swing.GroupLayout(panelPaths);
        panelPaths.setLayout(panelPathsLayout);
        panelPathsLayout.setHorizontalGroup(
            panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelPathsLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(panelPathsLayout.createSequentialGroup()
                        .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel12)
                            .addComponent(jLabel2))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(panelPathsLayout.createSequentialGroup()
                                .addComponent(fieldPheno, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(selPhenotypeBt))
                            .addGroup(panelPathsLayout.createSequentialGroup()
                                .addGap(2, 2, 2)
                                .addComponent(fieldGeno, javax.swing.GroupLayout.DEFAULT_SIZE, 452, Short.MAX_VALUE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(selGenotypeBt))))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, panelPathsLayout.createSequentialGroup()
                        .addComponent(labelOutputDir)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(selOutputDir)))
                .addContainerGap())
            .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(panelPathsLayout.createSequentialGroup()
                    .addGap(147, 147, 147)
                    .addComponent(fieldOutputDir, javax.swing.GroupLayout.DEFAULT_SIZE, 454, Short.MAX_VALUE)
                    .addGap(111, 111, 111)))
        );
        panelPathsLayout.setVerticalGroup(
            panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelPathsLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(selOutputDir, javax.swing.GroupLayout.PREFERRED_SIZE, 19, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(labelOutputDir))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel12)
                    .addComponent(fieldGeno, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(selGenotypeBt, javax.swing.GroupLayout.PREFERRED_SIZE, 19, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel2)
                    .addComponent(fieldPheno, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(selPhenotypeBt, javax.swing.GroupLayout.PREFERRED_SIZE, 19, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            .addGroup(panelPathsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                .addGroup(panelPathsLayout.createSequentialGroup()
                    .addGap(15, 15, 15)
                    .addComponent(fieldOutputDir, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addContainerGap(71, Short.MAX_VALUE)))
        );

        panelFilesTitle.add(panelPaths, java.awt.BorderLayout.CENTER);

        panelParametersTitle.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        panelParametersTitle.setLayout(new java.awt.BorderLayout());

        jLabel1.setBackground(javax.swing.UIManager.getDefaults().getColor("Button.disabledToolBarBorderBackground"));
        jLabel1.setText("GWAS Parameters:");
        jLabel1.setOpaque(true);
        panelParametersTitle.add(jLabel1, java.awt.BorderLayout.PAGE_START);

        panelParameters.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));

        jLabel8.setText("GWAS Model:");

        fieldModel.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "Full", "Naive" }));

        jLabel3.setText("Significance Level:");

        fieldSignificance.setText("0.05");

        jLabel4.setText("Correction Method:");

        fieldCorrection.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "Bonferroni", "FDR" }));

        jLabel5.setText("Number of Best SNPs:");

        fieldSNPs.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10" }));
        fieldSNPs.setSelectedIndex(9);

        jLabel6.setText("Filtering:");

        fieldFiltering.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "TRUE", "FALSE" }));
        fieldFiltering.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                fieldFilteringActionPerformed(evt);
            }
        });

        labelPloidy.setText("Ploidy:");

        fieldPloidy.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "4", "2" }));

        javax.swing.GroupLayout panelParametersLayout = new javax.swing.GroupLayout(panelParameters);
        panelParameters.setLayout(panelParametersLayout);
        panelParametersLayout.setHorizontalGroup(
            panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelParametersLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(panelParametersLayout.createSequentialGroup()
                        .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel8)
                            .addComponent(jLabel3)
                            .addComponent(labelPloidy))
                        .addGap(28, 28, 28)
                        .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(fieldModel, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(fieldSignificance)
                            .addComponent(fieldPloidy, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                    .addGroup(panelParametersLayout.createSequentialGroup()
                        .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, panelParametersLayout.createSequentialGroup()
                                .addComponent(jLabel5)
                                .addGap(5, 5, 5))
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, panelParametersLayout.createSequentialGroup()
                                .addComponent(jLabel6)
                                .addGap(95, 95, 95))
                            .addGroup(panelParametersLayout.createSequentialGroup()
                                .addComponent(jLabel4)
                                .addGap(22, 22, 22)))
                        .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(fieldCorrection, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(fieldSNPs, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(fieldFiltering, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))))
                .addContainerGap(74, Short.MAX_VALUE))
        );
        panelParametersLayout.setVerticalGroup(
            panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, panelParametersLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fieldPloidy, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(labelPloidy))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel8)
                    .addComponent(fieldModel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(6, 6, 6)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel3)
                    .addComponent(fieldSignificance, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(4, 4, 4)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel4)
                    .addComponent(fieldCorrection, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(6, 6, 6)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(fieldSNPs, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(panelParametersLayout.createSequentialGroup()
                        .addGap(4, 4, 4)
                        .addComponent(jLabel5)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(panelParametersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel6)
                    .addComponent(fieldFiltering, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(34, 34, 34))
        );

        panelParametersTitle.add(panelParameters, java.awt.BorderLayout.CENTER);

        panelFiltersTitle.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        panelFiltersTitle.setLayout(new java.awt.BorderLayout());

        panelFilters.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));

        fieldFilterMIND.setText("0.1");
        fieldFilterMIND.setPreferredSize(new java.awt.Dimension(40, 15));

        fieldFilterGENO.setText("0.1");
        fieldFilterGENO.setPreferredSize(new java.awt.Dimension(40, 15));
        fieldFilterGENO.setRequestFocusEnabled(false);

        fieldFilterHWE.setText("1e-10");
        fieldFilterHWE.setPreferredSize(new java.awt.Dimension(40, 15));
        fieldFilterHWE.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                fieldFilterHWEActionPerformed(evt);
            }
        });

        fieldFilterMAF.setText("0.01");
        fieldFilterMAF.setPreferredSize(new java.awt.Dimension(25, 15));
        fieldFilterMAF.setRequestFocusEnabled(false);

        jLabel7.setText("<html>Minimum Minor Allele Frequency (MAF) for a SNP to be kept: </html>");
        jLabel7.setName("<html>Exclude SNPs with minor <br> allele frequency less than </html>"); // NOI18N

        jLabel15.setText("<html>Maximum proportion of missing values for a SNP to be kept </html>");
        jLabel15.setName("<html>Exclude SNPs with minor <br> allele frequency less than </html>"); // NOI18N

        jLabel16.setText("<html>Maximum proportion of missing values for a sample to be kept: </html>");
        jLabel16.setName("<html>Exclude SNPs with minor <br> allele frequency less than </html>"); // NOI18N
        jLabel16.setPreferredSize(new java.awt.Dimension(440, 15));

        jLabel17.setText("<html>Filters out SNPs with HWE exact test p-value below threshold: </html>");
        jLabel17.setName("<html>Exclude SNPs with minor <br> allele frequency less than </html>"); // NOI18N

        javax.swing.GroupLayout panelFiltersLayout = new javax.swing.GroupLayout(panelFilters);
        panelFilters.setLayout(panelFiltersLayout);
        panelFiltersLayout.setHorizontalGroup(
            panelFiltersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelFiltersLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelFiltersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(panelFiltersLayout.createSequentialGroup()
                        .addComponent(jLabel16, javax.swing.GroupLayout.PREFERRED_SIZE, 232, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fieldFilterGENO, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(panelFiltersLayout.createSequentialGroup()
                        .addComponent(jLabel17, javax.swing.GroupLayout.PREFERRED_SIZE, 232, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fieldFilterHWE, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(panelFiltersLayout.createSequentialGroup()
                        .addComponent(jLabel15, javax.swing.GroupLayout.PREFERRED_SIZE, 232, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fieldFilterMIND, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(panelFiltersLayout.createSequentialGroup()
                        .addComponent(jLabel7, javax.swing.GroupLayout.PREFERRED_SIZE, 232, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fieldFilterMAF, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        panelFiltersLayout.setVerticalGroup(
            panelFiltersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelFiltersLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelFiltersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel7, javax.swing.GroupLayout.PREFERRED_SIZE, 40, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(fieldFilterMAF, javax.swing.GroupLayout.PREFERRED_SIZE, 40, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(panelFiltersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(jLabel15)
                    .addComponent(fieldFilterMIND, javax.swing.GroupLayout.PREFERRED_SIZE, 40, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(panelFiltersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(jLabel16, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(fieldFilterGENO, javax.swing.GroupLayout.PREFERRED_SIZE, 40, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(panelFiltersLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(fieldFilterHWE, javax.swing.GroupLayout.PREFERRED_SIZE, 40, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel17))
                .addGap(22, 22, 22))
        );

        jLabel7.getAccessibleContext().setAccessibleName("");

        panelFiltersTitle.add(panelFilters, java.awt.BorderLayout.CENTER);

        jLabel18.setBackground(javax.swing.UIManager.getDefaults().getColor("Button.disabledToolBarBorderBackground"));
        jLabel18.setText("Quality Control Filters:");
        jLabel18.setOpaque(true);
        panelFiltersTitle.add(jLabel18, java.awt.BorderLayout.PAGE_START);

        javax.swing.GroupLayout panelInputsLayout = new javax.swing.GroupLayout(panelInputs);
        panelInputs.setLayout(panelInputsLayout);
        panelInputsLayout.setHorizontalGroup(
            panelInputsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelInputsLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(panelInputsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                    .addComponent(panelFilesTitle, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(panelInputsLayout.createSequentialGroup()
                        .addComponent(panelParametersTitle, javax.swing.GroupLayout.PREFERRED_SIZE, 352, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(panelFiltersTitle, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(17, Short.MAX_VALUE))
        );
        panelInputsLayout.setVerticalGroup(
            panelInputsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(panelInputsLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(panelFilesTitle, javax.swing.GroupLayout.PREFERRED_SIZE, 121, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(panelInputsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(panelFiltersTitle, javax.swing.GroupLayout.PREFERRED_SIZE, 234, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(panelParametersTitle, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE))
                .addGap(218, 218, 218))
        );

        add(panelInputs, java.awt.BorderLayout.CENTER);
    }// </editor-fold>//GEN-END:initComponents

    private void selOutputDirActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_selOutputDirActionPerformed
        if (evt.getSource() == selOutputDir) {
            fc.setCurrentDirectory(getLastDir());
            fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
            int returnVal = fc.showOpenDialog(ViewInputs.this);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                File file = fc.getSelectedFile();
                prefs.put(LAST_USED_FOLDER, fc.getSelectedFile().getAbsolutePath());
                //This is where a real application would open the file.
                fieldOutputDir.setText(file.getAbsolutePath());
                setEnabledInputs(true);
            }
        }
    }//GEN-LAST:event_selOutputDirActionPerformed

    private void selPhenotypeBtActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_selPhenotypeBtActionPerformed
        if (evt.getSource() == selPhenotypeBt) {
            fc.setCurrentDirectory(getLastDir());
            int returnVal = fc.showOpenDialog(ViewInputs.this);
            fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            if (returnVal == JFileChooser.APPROVE_OPTION) {
                File file = fc.getSelectedFile();
                //This is where a real application would open the file.
                fieldPheno.setText(file.getAbsolutePath());
            }
        }
    }//GEN-LAST:event_selPhenotypeBtActionPerformed

    private void selGenotypeBtActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_selGenotypeBtActionPerformed
        //Handle open button action.
        if (evt.getSource() == selGenotypeBt) {
            fc.setCurrentDirectory(getLastDir());
            fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
            int returnVal = fc.showOpenDialog(ViewInputs.this);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                File file = fc.getSelectedFile();
                //This is where a real application would open the file.
                fieldGeno.setText(file.getAbsolutePath());
                prefs.put(LAST_USED_FOLDER, fc.getSelectedFile().getParent());
            }
        }        // TODO add your handling code here:
    }//GEN-LAST:event_selGenotypeBtActionPerformed

    private void fieldFilteringActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_fieldFilteringActionPerformed
        boolean flag = false;
        System.out.println(fieldFiltering.getSelectedItem().toString());
        if (fieldFiltering.getSelectedItem().toString().equals("FALSE")) {
            flag = false;
        } else {
            flag = true;
        }

        for (Component c : panelFilters.getComponents())
            c.setEnabled(flag);
    }//GEN-LAST:event_fieldFilteringActionPerformed

    private void fieldFilterHWEActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_fieldFilterHWEActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_fieldFilterHWEActionPerformed


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JComboBox<String> fieldCorrection;
    private javax.swing.JTextField fieldFilterGENO;
    private javax.swing.JTextField fieldFilterHWE;
    private javax.swing.JTextField fieldFilterMAF;
    private javax.swing.JTextField fieldFilterMIND;
    private javax.swing.JComboBox<String> fieldFiltering;
    private javax.swing.JTextField fieldGeno;
    private javax.swing.JComboBox<String> fieldModel;
    private javax.swing.JTextField fieldOutputDir;
    private javax.swing.JTextField fieldPheno;
    private javax.swing.JComboBox<String> fieldPloidy;
    private javax.swing.JComboBox<String> fieldSNPs;
    private javax.swing.JTextField fieldSignificance;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel12;
    private javax.swing.JLabel jLabel14;
    private javax.swing.JLabel jLabel15;
    private javax.swing.JLabel jLabel16;
    private javax.swing.JLabel jLabel17;
    private javax.swing.JLabel jLabel18;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JLabel jLabel8;
    private javax.swing.JLabel labelOutputDir;
    private javax.swing.JLabel labelPloidy;
    private javax.swing.JPanel panelFilesTitle;
    private javax.swing.JPanel panelFilters;
    private javax.swing.JPanel panelFiltersTitle;
    private javax.swing.JPanel panelInputs;
    private javax.swing.JPanel panelParameters;
    private javax.swing.JPanel panelParametersTitle;
    private javax.swing.JPanel panelPaths;
    private javax.swing.JButton selGenotypeBt;
    private javax.swing.JButton selOutputDir;
    private javax.swing.JButton selPhenotypeBt;
    // End of variables declaration//GEN-END:variables
}
