# coding: utf-8

""" Matplotlib figure functionality for the Traits GUI. """

from __future__ import division, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Third-party
import wx
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.pyplot import Figure, subplot
from matplotlib.ticker import MaxNLocator

# Editor factories for matplotlib figures
from traitsui.wx.editor import Editor
from traitsui.wx.basic_editor_factory import BasicEditorFactory


class _MPLFigureEditor(Editor):
    """ Editor class for containing a matplotlib figure within a Traits GUI. """

    scrollable  = True

    def init(self, parent):
        self.control = self._create_canvas(parent)
        self.set_tooltip()

    def update_editor(self):
        pass

    def _create_canvas(self, parent):
        """ Create the matplotlib canvas """

        panel = wx.Panel(parent, -1, style=wx.CLIP_CHILDREN)
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(sizer)
        
        mpl_control = FigureCanvasWxAgg(panel, -1, self.value)
        sizer.Add(mpl_control, 1, wx.LEFT | wx.TOP | wx.GROW)
        toolbar = NavigationToolbar2Wx(mpl_control)
        sizer.Add(toolbar, 0, wx.EXPAND)
        self.value.canvas.SetMinSize((10,10))        
        return panel


class MPLFigureEditor(BasicEditorFactory):
    """ Factory class for generating editors that contain matplotlib figures and
        can be placed within a Traits GUI. """
    klass = _MPLFigureEditor