# style.py

def get_stylesheet():
    """
    Returns the QSS stylesheet string using a non-default font.
    :return: QSS string.
    """
    stylesheet = """
        * {
            font-family: 'Roboto';  /* Use any non-default font, e.g., Helvetica or Arial */
        }

        QMainWindow {
            background-color: #2b2b2b;
        }

        QFrame {
            background-color: #1c1c1c;
            border: 2px solid #3c3c3c;
            border-radius: 10px;
            min-width: 300px;  /* Set a minimum width */
            min-height: 200px; /* Set a minimum height */
        }
        QPushButton {
            background-color: #3c3c3c;
            color: #ffffff;
        }

        QPushButton:hover {
            background-color: #444444;
            border-color: #6a6a6a;
        }

        QLabel {
            color: #ffffff;
            font-size: 16px;
        }

        QMenuBar {
            background-color: #2b2b2b;
            color: #ffffff;
        }

        QMenuBar::item:selected {
            background: #444444;
        }

        QMenu {
            background-color: #2b2b2b;
            color: #ffffff;
        }

        QMenu::item:selected {
            background-color: #444444;
        }

        QStatusBar {
            background-color: #2b2b2b;
            color: #ffffff;
        }

        QComboBox {
            background-color: #3c3c3c;
            color: #ffffff;
            border: 1px solid #5c5c5c;
            padding: 5px;
        }

        QLineEdit {
            background-color: #2b2b2b;
            color: #ffffff;
            border: 1px solid #5c5c5c;
            padding: 5px;
        }

        QProgressBar {
            background-color: #2b2b2b;
            border: 1px solid #5c5c5c;
            color: #ffffff;
        }

        QProgressBar::chunk {
            background-color: #444444;
            width: 20px;
        }
    """
    return stylesheet

# dark_theme.py

def apply_dark_theme(axes):
    """
    Apply dark theme styling to a Matplotlib axes object.
    
    :param axes: Matplotlib axes object to style
    """
    # Set background color for the figure and axes
    axes.figure.patch.set_facecolor('#1a1a2e')  # Dark blue background for figure
    axes.set_facecolor('#2b2b2b')  # Slightly lighter background for the plot area

    # Set axes spines color (borders around the plot)
    axes.spines['bottom'].set_color('#ffffff')  # White for bottom axis
    axes.spines['top'].set_color('#ffffff')  # White for top axis
    axes.spines['left'].set_color('#ffffff')  # White for left axis
    axes.spines['right'].set_color('#ffffff')  # White for right axis

    # Set axis labels and tick colors
    axes.xaxis.label.set_color('#ffffff')  # Set X axis label color to white
    axes.yaxis.label.set_color('#ffffff')  # Set Y axis label color to white
    axes.tick_params(axis='x', colors='#ffffff')  # Set X axis tick color to white
    axes.tick_params(axis='y', colors='#ffffff')  # Set Y axis tick color to white

    # Set gridlines color
    axes.grid(True, color='#444444')  # Light gray grid lines

    # Set title and text color
    axes.title.set_color('#ffffff')  # Set plot title color to white
