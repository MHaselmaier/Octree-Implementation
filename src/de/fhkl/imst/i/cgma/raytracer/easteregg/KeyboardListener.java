package de.fhkl.imst.i.cgma.raytracer.easteregg;

import java.awt.event.KeyEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.util.Arrays;


public class KeyboardListener implements java.awt.event.KeyListener
{
    private final int[] last_ten_pressed_keys = new int[10];
    private int counter = 0;
    private Pattern pattern;
    private PropertyChangeSupport propertyChangeSupport;

    public KeyboardListener(Pattern pattern)
    {
        this.propertyChangeSupport = new PropertyChangeSupport(this);
        this.pattern = pattern;
    }

    public void addPropertyChangeListener(PropertyChangeListener propertyChangeListener)
    {
        this.propertyChangeSupport.addPropertyChangeListener(propertyChangeListener);
    }

    public void removePropertyChangeListener(PropertyChangeListener propertyChangeListener)
    {
        this.propertyChangeSupport.removePropertyChangeListener(propertyChangeListener);
    }

    @Override
    public void keyPressed(KeyEvent keyEvent)
    {
        this.last_ten_pressed_keys[this.counter++] = keyEvent.getKeyCode();
        if(10 == this.counter)
        {
            counter = 0;
            this.propertyChangeSupport.firePropertyChange("match_pattern", false, checkPressedKeyPattern());
        }
    }

    private boolean checkPressedKeyPattern()
    {
        return Arrays.equals(this.last_ten_pressed_keys, this.pattern.getPattern());
    }

    @Override
    public void keyTyped(KeyEvent e)
    {
    }

    @Override
    public void keyReleased(KeyEvent e)
    {
    }
}
