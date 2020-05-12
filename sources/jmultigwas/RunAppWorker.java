/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jmultigwas;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import javax.swing.SwingWorker;

/**
 *
 * @author lg
 */
public class RunAppWorker extends SwingWorker<Void, String> {

    Controller controller;
    String configFilenamePath;
    String outputDir;

    public RunAppWorker(String configFilenamePath, String outputDir, Controller controller) {
        this.controller = controller;
        this.configFilenamePath = configFilenamePath;
        this.outputDir = outputDir;
        
    }

    @Override
    protected Void doInBackground() throws Exception {
        ProcessBuilder processBuilder = new ProcessBuilder();
        processBuilder.directory(new File(outputDir));
        
        // -- Linux ... Run a shell command

        //processBuilder.command("bash", "-c", "pwd");
        //System.out.println ("....Out");
        //System.exit (0);
        Path outputPath       = Paths.get(outputDir);
        String outputDirName  = outputPath.getFileName().toString();
        String configFilename = outputDirName + ".config";
        
        String commandString = "multiGWAS " + configFilename;
        
        processBuilder.command("bash", "-c", commandString, outputPath.toString());
        System.out.println (">>> Command: " + processBuilder.command());
        try {
            Process process = processBuilder.start();
            InputStreamReader isr = new InputStreamReader(process.getInputStream());
            int c;
            while ((c = isr.read()) >= 0) {
                controller.sendToOutputs((char) c);
                System.out.print((char) c);
                //System.out.flush();
            }
            System.out.println("End of execution");
            int exitVal = process.waitFor();
            if (exitVal == 0) {
                System.out.println("Success!");
            } else {
                System.out.println("Abnormal!");
                System.exit(1);
            }
            controller.onEndOfExecution();

        } catch (IOException e) {
            e.printStackTrace();
            //} catch (InterruptedException e) {
            //  e.printStackTrace();
        }
        return null;
    }

    @Override
    protected void process(List<String> chunks) {
        
    }

}
