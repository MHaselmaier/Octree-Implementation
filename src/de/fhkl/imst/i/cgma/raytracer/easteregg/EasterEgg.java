package de.fhkl.imst.i.cgma.raytracer.easteregg;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

public abstract class EasterEgg implements PropertyChangeListener
{
    private final KeyboardListener keyboardListener;
    public abstract void trigger();

    public EasterEgg(Pattern pattern)
    {
        this.keyboardListener = new KeyboardListener(pattern);
        this.keyboardListener.addPropertyChangeListener(this);
    }

    public KeyboardListener getKeyboardListener()
    {
        return this.keyboardListener;
    }

    @Override
    public void propertyChange(PropertyChangeEvent propertyChangeEvent)
    {
        this.trigger();
    }
}
