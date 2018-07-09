package de.fhkl.imst.i.cgma.raytracer.easteregg;

import javax.swing.*;

public class KonamiCodeDemo
{
    public static void main(String[] args)
    {
        JFrame mainFrame = new JFrame();
        mainFrame.setTitle("Demo");
        mainFrame.setSize(400,400);
        mainFrame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        mainFrame.addKeyListener(new AllYourBaseAreBelongToUsEasterEgg(new KonamiCode()).getKeyboardListener());
        mainFrame.setVisible(true);
    }

}