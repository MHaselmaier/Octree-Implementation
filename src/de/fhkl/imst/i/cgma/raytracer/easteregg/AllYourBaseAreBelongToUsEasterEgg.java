package de.fhkl.imst.i.cgma.raytracer.easteregg;

import javafx.application.Platform;
import javafx.embed.swing.JFXPanel;
import javafx.scene.Scene;
import javafx.scene.layout.StackPane;
import javafx.scene.media.Media;
import javafx.scene.media.MediaPlayer;
import javafx.scene.media.MediaView;

import javax.swing.*;
import java.io.File;
import java.net.MalformedURLException;
import java.net.URL;


public class AllYourBaseAreBelongToUsEasterEgg extends  EasterEgg
{
    private static final File VIDEO_FILE = new File("res/AYBABTU.mp4");
    private URL mediaFileURL;

    public AllYourBaseAreBelongToUsEasterEgg(Pattern pattern)
    {
        super(pattern);
        try
        {
            this.mediaFileURL = AllYourBaseAreBelongToUsEasterEgg.VIDEO_FILE.toURI().toURL();
        }
        catch (MalformedURLException malFormedURLException)
        {
            malFormedURLException.printStackTrace();
        }
    }

    @Override
    public void trigger()
    {
        final JFrame mediaPlayerFrame = createJFrame();
        final JFXPanel javaFXPanel = createJFXPanel();
        mediaPlayerFrame.add(javaFXPanel);
        mediaPlayerFrame.setVisible(true);
        Platform.runLater(() -> initFX(mediaPlayerFrame, javaFXPanel));
    }

    private JFrame createJFrame()
    {
        JFrame jFrame = new JFrame();
        jFrame.setSize(400,400);
        return jFrame;
    }

    private JFXPanel createJFXPanel()
    {
        JFXPanel jfxPanel = new JFXPanel();
        jfxPanel.setSize(400,400);
        return jfxPanel;
    }

    private void initFX(JFrame jFrame, JFXPanel jfxPanel)
    {
        StackPane root = new StackPane();
        root.setMaxSize(400,400);
        jfxPanel.setScene(new Scene(root, 400,400));
        MediaPlayer mediaPlayer = createMediaPlayer(jFrame);
        MediaView mediaView = createMediaView(mediaPlayer, 400,400);
        root.getChildren().add(mediaView);
        jfxPanel.setVisible(true);
        mediaPlayer.play();
    }

    private MediaPlayer createMediaPlayer(JFrame jFrame)
    {
        MediaPlayer mediaPlayer = new MediaPlayer(new Media(this.mediaFileURL.toString()));
        mediaPlayer.setOnEndOfMedia(() -> jFrame.setVisible(false));
        return mediaPlayer;
    }

    private MediaView createMediaView(MediaPlayer mediaPlayer, int height, int width)
    {
        MediaView mediaView = new MediaView(mediaPlayer);
        mediaView.setFitHeight(height);
        mediaView.setFitWidth(width);
        return mediaView;
    }
}