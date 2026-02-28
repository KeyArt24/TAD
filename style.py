


scroll_vertical = """
    QScrollBar:vertical
    {
        background-color: white;
        width: 15px;
        margin: 15px 3px 15px 3px;
        border: 1px transparent #2A2929;
        border-radius: 4px;
    }

    QScrollBar::handle:vertical
    {
        background-color: gray;         /* #605F5F; */
        min-height: 5px;
        border-radius: 4px;
    }

    QScrollBar::sub-line:vertical
    {
        margin: 3px 0px 3px 0px;
        border-image: url(:/qss_icons/rc/up_arrow_disabled.png);
        height: 10px;
        width: 10px;
        subcontrol-position: top;

    }

    QScrollBar::add-line:vertical
    {
        margin: 3px 0px 3px 0px;
        border-image: url(:/qss_icons/rc/down_arrow_disabled.png);
        height: 10px;
        width: 10px;
        subcontrol-position: bottom;
        subcontrol-origin: margin;
    }

    QScrollBar::up-arrow:vertical, QScrollBar::down-arrow:vertical
    {
        background: none;
    }

    QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical
    {
        background: none;
    }
"""

scroll_horizontal = """
    QScrollBar:horizontal {
        border: 1px transparent #2A2929;
        background: gray;
        height: 10px;
        border-radius: 4px;

    }
    
    QScrollBar::handle:horizontal
    {
        background-color: lightgray;         /* #605F5F; */
        min-width: 5px;
        border-radius: 4px;
    }
    
    QScrollBar::sub-line:width
    {
        margin: 3px 0px 3px 0px;
        border-image: url(:/qss_icons/rc/up_arrow_disabled.png);
        height: 10px;
        width: 10px;
        subcontrol-position: top;

    }
    
    QScrollBar::add-line:horizontal
    {
        margin: 3px 0px 3px 0px;
        border-image: url(:/qss_icons/rc/down_arrow_disabled.png);
        height: 10px;
        width: 10px;
        subcontrol-position: bottom;
        subcontrol-origin: margin;
    }
    
    QScrollBar::up-arrow:horizontal, QScrollBar::down-arrow:horizontal
    {
        background: none;
    }

    QScrollBar::add-page:horizontal, QScrollBar::sub-page:horizontal
    {
        background: none;
    }
    
    """

tab_lay = '''
            QTabBar::tab-bar {alignment: center}
            QTabBar::tab  {background-color:rgb(220, 220, 220); border-radius: 4px; height: 20px; width: 150px; margin: 2px;}
            QTabWidget::pane {position; absolute}
            '''

css_buttons_sample = "background-color: white; color: black"