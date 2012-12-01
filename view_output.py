import wx
import  wx.lib.mixins.listctrl as listmix
import sys
import csv


class MyListCtrl(listmix.ListCtrlAutoWidthMixin, wx.ListCtrl):
    def __init__(self, p, my_id, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        wx.ListCtrl.__init__(self, p, my_id,
                             pos, size, style)
        listmix.ListCtrlAutoWidthMixin.__init__(self)


class ViewOutputFrame(wx.Frame, listmix.ColumnSorterMixin):
    def __init__(self, parent, id, string, output_file_to_load):
        # wx.Dialog.__init__(self, parent, -1, style=wx.WANTS_CHARS)
        wx.Frame.__init__(self, parent, -1, string)
        # ToolBar at the top of the window
        toolbar = wx.ToolBar(self, -1, style=wx.TB_HORIZONTAL | wx.NO_BORDER)
        toolbar.SetToolBitmapSize(size=(24, 24))
        toolbar.AddSimpleTool(wx.ID_OPEN, self.get_bmp(wx.ART_FILE_OPEN),
            "Load", "Load a text file")
        toolbar.AddSimpleTool(wx.ID_SAVE, self.get_bmp(wx.ART_FILE_SAVE),
            "Save", "Save the text file")
        toolbar.AddSimpleTool(wx.ID_ABOUT, self.get_bmp(wx.ART_INFORMATION),
            "About", "About message")
        toolbar.AddSeparator()
        toolbar.AddSimpleTool(wx.ID_EXIT, self.get_bmp(wx.ART_QUIT),
            "Exit", "Exit the program")
        toolbar.Realize()
        self.SetToolBar(toolbar)

        tID = wx.NewId()

        sizer = wx.BoxSizer(wx.VERTICAL)

        self.list = MyListCtrl(self, tID,
                               style=wx.LC_REPORT
                               | wx.LC_SORT_ASCENDING
                               | wx.LC_EDIT_LABELS
                               | wx.BORDER_NONE)

        # self.list.SetImageList(self.il, wx.IMAGE_LIST_SMALL)
        sizer.Add(self.list, 1, wx.EXPAND)

        self.read_output_file(output_file_to_load)

        # Now that the list exists we can init the other base class,
        # see wx/lib/mixins/listctrl.py
        listmix.ColumnSorterMixin.__init__(self, 3)
        #self.SortListItems(0, True)

        self.SetSizer(sizer)
        self.SetAutoLayout(True)

    def read_output_file(self, output_file):
        # define the columns
        self.list.InsertColumn(0, 'ID')
        self.list.InsertColumn(1, 'Target Coordinate', wx.LIST_FORMAT_RIGHT)
        self.list.InsertColumn(2, 'Primer Coordinates')
        self.list.InsertColumn(3, 'Psi Target')
        self.list.InsertColumn(4, 'Upstream Primer')
        self.list.InsertColumn(5, 'Downstream Primer')
        self.list.InsertColumn(6, 'Average TM')
        self.list.InsertColumn(7, 'Skipping Prod. Size')
        self.list.InsertColumn(8, 'Inc. Prod. Sized')
        self.list.InsertColumn(9, 'Upstream Exon Coord.')
        self.list.InsertColumn(10, 'Upstream Psi')
        self.list.InsertColumn(11, 'Downstream Exon Coord.')
        self.list.InsertColumn(12, 'Downstream Psi')

        # get information from results file
        with open(output_file) as handle:
            self.results = list(csv.reader(handle, delimiter='\t'))[1:]

        for i, data in enumerate(self.results):
            # index = self.list.InsertImageStringItem(sys.maxint, data[0], self.idx1)
            index = self.list.InsertStringItem(sys.maxint, data[0])
            for col in range(1, len(data)):
                self.list.SetStringItem(index, col, data[col])
            # self.list.SetStringItem(index, 1, data[1])
            # self.list.SetStringItem(index, 2, data[2])
            self.list.SetItemData(index, i)

        self.list.SetColumnWidth(0, wx.LIST_AUTOSIZE)
        self.list.SetColumnWidth(1, wx.LIST_AUTOSIZE)
        self.list.SetColumnWidth(2, 100)

        # show how to select an item
        self.list.SetItemState(5, wx.LIST_STATE_SELECTED, wx.LIST_STATE_SELECTED)

        self.current_selection = 0

    def get_bmp(self, pic_id):
        return wx.ArtProvider.GetBitmap(pic_id, wx.ART_TOOLBAR, wx.Size(24, 24))

    def GetListCtrl(self):
        return self.list

    # Used by the ColumnSorterMixin, see wx/lib/mixins/listctrl.py

    def GetSortImages(self):
        return (self.sm_dn, self.sm_up)


class TestViewOutputApp(wx.App):
    '''
    This is just a test app to see if this works indepently
    '''
    def OnInit(self):
        wx.InitAllImageHandlers()
        view_output = ViewOutputFrame(None, -1, "", 'output.txt')
        view_output.Show()
        return 1

if __name__ == '__main__':
    ViewOutputApp = TestViewOutputApp(0)
    ViewOutputApp.MainLoop()


