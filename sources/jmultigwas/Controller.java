package jmultigwas;

import java.awt.BorderLayout;
import java.awt.Desktop;
import java.awt.Dimension;
import java.io.File;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Paths;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
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

        tabInputs = new ViewInputs(this, model);
        tabOutputs = new ViewOutputs(size);
        tabResults = new ViewResults(size);
        tabFiles = new ViewFiles(this);
        tabOutputs.init();

        viewTabs.addTab("Inputs", tabInputs);
        viewTabs.addTab("Outputs", tabOutputs);
        viewTabs.addTab("Results", tabResults);
        viewTabs.addTab("Files", tabFiles);

        //onNewProject();
    }

    public String getToolsToRun() {
        String tools = viewToolBar.getToolsToRun();
        

        return (tools);
    }

    public void onDefaultButton() {
        tabInputs.setDefaults();
    }

    public void onRunApplication() {
        if (tabInputs.checkCompleteInfo()==false) {
            JOptionPane.showMessageDialog(this, "Incomplete information", "MultiGWAS warning",
                    JOptionPane.WARNING_MESSAGE);
        }else {
            String outputDir  = tabInputs.getOutputDir();
            String values     = tabInputs.getInputValues();

            model.runApplication(outputDir, values);
            viewTabs.setSelectedIndex(1);
        }
    }

    public void onCancelButton() {
        viewTabs.setSelectedComponent(tabInputs);
    }

    public void onEndOfExecution() {
        String workingDir = tabInputs.getOutputDir();
        String dirName = Paths.get(workingDir).getFileName().toString();
        String outputDir = workingDir + File.separator + "out-" + dirName;
        String reportDir = outputDir + File.separator + "report";
        String htmlFilename = outputDir + File.separator + "multiGWAS-report.html";

        browseFile(htmlFilename, reportDir);
    }

    public void browseFile(String url, String reportDir) {
        String myOS = System.getProperty("os.name").toLowerCase();
        OUT("(Your operating system is: " + myOS + ")\n");

        try {
            if (myOS.contains("win") && Desktop.isDesktopSupported()) { // Probably Windows
                OUT(" -- Going with Desktop.browse ...");
                Desktop desktop = Desktop.getDesktop();
                desktop.browse(new URI(url));
            } else { // Definitely Non-windows
                Runtime runtime = Runtime.getRuntime();
                if (myOS.contains("mac")) { // Apples
                    OUT(" -- Going on Apple with 'open'...");
                    runtime.exec("open " + url);
                } else if (myOS.contains("nix") || myOS.contains("nux")) { // Linux flavours 
                    OUT(" -- Going on Linux with 'xdg-open'...");
                    runtime.exec("xdg-open " + url);
                } else {
                    OUT("I was unable/unwilling to launch a browser in your OS :( #SadFace");
                }
            }
            OUT("\nThings have finished.\nI hope you're OK.");
        } catch (IOException | URISyntaxException eek) {
            OUT("**Stuff wrongly: " + eek.getMessage());
        }

        tabResults.showResults(url);
        viewTabs.setSelectedIndex(2);
        OUT(">>> Report file: " + url);
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

    public void sendToOutputs(String s) {
        int n = s.length();
        for (int i = 0; i < n; i++) {
            tabOutputs.append(s.charAt(i));
        }
    }

    private void OUT(String string) {
        System.out.println(string);
    }

}
