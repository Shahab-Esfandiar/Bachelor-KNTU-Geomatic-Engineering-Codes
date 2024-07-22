import arcpy
import pythonaddins

def MyResetButton(toolID,defaultValue):
    toolID.value = defaultValue

class ButtonClass1(object):
    """Implementation for MyAddin_addin.button (Button)"""
    def __init__(self):
        self.enabled = True
        self.checked = False
    def onClick(self):
        arcpy.SelectLayerByAttribute_management(mylayer,"CLEAR_SELECTION")
        self.MyResetButton(MyAddin_addin.combobox,"")
        MyAddin_addin.combobox.enabled = False
        arcpy.RefreshActiveView()

class ComboBoxClass2(object):
    """Implementation for MyAddin_addin.combobox (ComboBox)"""
    def __init__(self):
        self.editable = True
        self.enabled = True
        self.dropdownWidth = 'WWWWWW'
        self.width = 'WWWWWW'
    def onSelChange(self, selection):
        arcpy.SelectLayerByAttribute_management(mylayer,"NEW_SELECTION","LAYER = '"+ selection + "'")
        arcpy.RefreshActiveView()
        mycount = int(arcpy.GetCount_management(mylayer).getOutput(0))
        pythonaddins.MessageBox(mycount,'The number of selected values', 0)
        
    def onEditChange(self, text):
        pass
    def onFocus(self, focused):
        global mylayer
        self.mxd=arcpy.mapping.MapDocument('current')
        mylayer=arcpy.mapping.ListLayers(self.mxd)[0]
        self.items=[]
        with arcpy.da.SearchCursor(mylayer, ["LAYER"]) as d:
            for row in d:
                self.items.append(row[0])
    def onEnter(self):
        pass
    def refresh(self):
        pass
