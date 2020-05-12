package jmultigwas;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.io.File;
import java.nio.file.Paths;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JTabbedPane;

class Controller extends JFrame {

    // Attributes
    Model model;
    JTabbedPane viewTabs;
    ViewInputs tabInputs;
    ViewFiles tabFiles;
    ViewToolBar viewToolBar;

    ViewOutputs tabOutputs;
    ViewResults tabResults;
    //JPanel tabInputs;

    JMenu menu, submenu;
    JMenuItem itemNew, itemOpen;

    // Methods
    public Controller(String text) {
        super(text);
        model = new Model(this);
        viewToolBar = new ViewToolBar(this);
        setTitle("JMultiGWAS Tool for GWAS");
    }
    
    public void init() {
        this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        this.setSize(900, 500);

        initViewTabs();

        //this.setContentPane(viewTabs);
        this.add(viewTabs, BorderLayout.CENTER);
        this.add(viewToolBar, BorderLayout.WEST);

        // this.setMenu();
        this.setVisible(true);
        tabInputs.setEnabledInputs(false);
    }

    public void initViewTabs() {
        viewTabs = new JTabbedPane();
        Dimension size = this.getSize();

        tabInputs  = new ViewInputs (this, model);
        tabOutputs = new ViewOutputs(size);
        tabResults = new ViewResults(size);
        tabFiles   = new ViewFiles(this);
        tabOutputs.init();

        viewTabs.addTab("Inputs", tabInputs);
        viewTabs.addTab("Outputs", tabOutputs);
        viewTabs.addTab("Results", tabResults);
        viewTabs.addTab("Files", tabFiles);

        //onNewProject();
    }
    
    public void onDefaultButton () {
        tabInputs.setDefaults ();
    }
    public void onRunApplication () {
        viewTabs.setSelectedIndex(1);
    }
    
    public void onCancelButton () {
        viewTabs.setSelectedComponent(tabInputs);
    }

    public void onEndOfExecution() {
        String workingDir   = tabInputs.getOutputDir();
        String dirName      = Paths.get (workingDir).getFileName ().toString ();
        String outputDir    = workingDir + File.separator + "out-" + dirName;
        String reportDir    = outputDir + File.separator + "report";
        String htmlFilename = outputDir + File.separator + "multiGWAS-report.html";
       
        tabResults.showResults(htmlFilename);
        viewTabs.setSelectedIndex(2);
        System.out.println (">>> report dir: " + reportDir);
        tabFiles.changeDir(reportDir);
        //tabResults = new ViewResults("/tmp/multiGWAS-report.html");
    }

    public void onNewButton() {
        tabInputs.clearInputs();
        //tabInputs.add(tabInputs, BorderLayout.CENTER);
        //paintComponents(getGraphics());
    }

    public void sendToOutputs(char c) {
        tabOutputs.append(c);
    }

}
