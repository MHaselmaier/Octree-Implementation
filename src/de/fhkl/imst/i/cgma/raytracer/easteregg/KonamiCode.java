package de.fhkl.imst.i.cgma.raytracer.easteregg;

public class KonamiCode implements Pattern
{
            private static final int    UP = 38;
            private static final int    DOWN = 40;
            private static final int    LEFT = 37;
            private static final int    RIGHT = 39;
            private static final int    B = 66;
            private static final int    A = 65;
            private static final int[]   PATTERN =
            {
                KonamiCode.UP,
                KonamiCode.UP,
                KonamiCode.DOWN,
                KonamiCode.DOWN,
                KonamiCode.LEFT,
                KonamiCode.RIGHT,
                KonamiCode.LEFT,
                KonamiCode.RIGHT,
                KonamiCode.B,
                KonamiCode.A
            };

    @Override
    public int[] getPattern()
    {
        return KonamiCode.PATTERN;
    }
}
