package jmultigwas;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

public class ViewOutputs extends JPanel {
    Dimension dimension;
    JScrollPane scroll;
    
    public ViewOutputs (Dimension d) {
        super ();
        setLayout (new BorderLayout ());
        dimension = new Dimension (d.width, (int) (0.9*d.height));       
    }

    JTextArea outputArea = new JTextArea();

    public void init() {
        setBackground (Color.GRAY);
        outputArea.setEditable(false);
        
        scroll = new JScrollPane(outputArea,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        scroll.setPreferredSize(dimension);
        //System.exit (0);
        add(scroll);
        //this.outputArea = outputArea;
    }

    public void append(char c) {
        outputArea.append(c + "");
        outputArea.setCaretPosition(outputArea.getDocument().getLength());
    }
}
