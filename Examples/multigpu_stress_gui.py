"""Launch the stress colony with graphical configuration and live diagnostics.

Run from the checkout:
    python Examples/multigpu_stress_gui.py

Choose population settings, then select GPUs in CellModeller's device dialog.
Click Run in the main window to start. This is a standalone GUI launcher, not
a model for Load Model; the underlying model is multigpu_stress.py.
Requires the same PyQt5, OpenGL and OpenCL dependencies as CellModellerGUI.py.
"""

import os
from pathlib import Path
import sys


def configure(parent=None):
    from PyQt5.QtWidgets import (QDialog, QFormLayout, QLabel, QComboBox,
                                 QSpinBox, QDoubleSpinBox, QDialogButtonBox)
    dialog = QDialog(parent)
    dialog.setWindowTitle('Multi-GPU stress test — colony settings')
    layout = QFormLayout(dialog)
    intro = QLabel('Grow a dense colony to the selected capacity, then stop growth.\n'
                   'GPU selection and workload weights follow in the device dialog.')
    intro.setWordWrap(True)
    layout.addRow(intro)
    mode = QComboBox()
    mode.addItems(['Explicit cell target', 'Estimate from available memory'])
    layout.addRow('Capacity selection', mode)
    fields = {}

    def integer(key, label, default, low, high):
        widget = QSpinBox()
        widget.setRange(low, high)
        widget.setValue(default)
        fields[key] = widget
        layout.addRow(label, widget)
        return widget

    target = integer('MAX_CELLS', 'Target cells', 100000, 1, 8000000)
    initial = integer('INITIAL_CELLS', 'Initial cells', 256, 1, 100000)
    contacts = integer('MAX_CONTACTS', 'Contacts per cell', 32, 8, 1024)
    integer('SPECIES', 'Intracellular species', 4, 1, 256)
    integer('SEED', 'Random seed', 12345, 0, 2147483647)
    integer('REPORT_EVERY', 'Report every N steps', 10, 1, 10000)
    grid = integer('MAX_SQS', 'Spatial grid bins (0 = automatic)', 0, 0, 2147483647)
    grid.setSpecialValueText('Automatic')
    for key, label, default, low, high in (
        ('GROWTH_RATE', 'Growth rate', 1.0, 0.01, 5.0),
        ('MEMORY_FRACTION', 'Memory fraction for estimate', 0.5, 0.05, 0.8),
    ):
        widget = QDoubleSpinBox()
        widget.setRange(low, high)
        widget.setSingleStep(0.05)
        widget.setValue(default)
        fields[key] = widget
        layout.addRow(label, widget)

    note = QLabel('Memory estimates are targets, not measured capacity limits. '
                  'Host RAM and transfers can limit performance before GPUs fill. '
                  'Large initial populations can take a long time to initialize. '
                  'Pickle saving is disabled for this test.')
    note.setWordWrap(True)
    layout.addRow(note)
    error = QLabel()
    error.setWordWrap(True)
    layout.addRow(error)
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    layout.addRow(buttons)

    def capacity_mode():
        target.setEnabled(mode.currentIndex() == 0)
        fields['MEMORY_FRACTION'].setEnabled(mode.currentIndex() == 1)

    def accept():
        if mode.currentIndex() == 0:
            if initial.value() > target.value():
                error.setText('Initial cells must not exceed the target.')
                return
            if target.value() * contacts.value() * 8 >= 2**31:
                error.setText('Reduce the target or contacts to fit supported index limits.')
                return
        dialog.accept()

    mode.currentIndexChanged.connect(capacity_mode)
    capacity_mode()
    buttons.accepted.connect(accept)
    buttons.rejected.connect(dialog.reject)
    if dialog.exec_() != QDialog.Accepted:
        return None
    values = {key: str(widget.value()) for key, widget in fields.items()}
    if mode.currentIndex() == 1:
        values['MAX_CELLS'] = 'auto'
    if grid.value() == 0:
        del values['MAX_SQS']
    return values


def main():
    # Resolve relative to this script, even when launched outside the checkout.
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root))
    from PyQt5.QtWidgets import QApplication, QDialog, QVBoxLayout, QLabel, QMessageBox
    from PyQt5.QtCore import QTimer, Qt
    from PyQt5 import uic

    app = QApplication(sys.argv)
    settings = configure()
    if settings is None:
        return 0
    # The model reads these during setup and reset. Restore the original
    # environment when this launcher exits, including an old manual grid limit.
    keys = {'CM_STRESS_' + key for key in settings} | {'CM_STRESS_MAX_SQS'}
    previous = {key: os.environ.get(key) for key in keys}
    try:
        os.environ.pop('CM_STRESS_MAX_SQS', None)
        os.environ.update({'CM_STRESS_' + key: value for key, value in settings.items()})
        # Match the existing GUI launcher, using the checkout's UI definition.
        import CellModeller.GUI.Renderers
        from CellModeller.GUI.PyGLCMViewer import PyGLCMViewer  # Qt custom widget
        ui = uic.loadUi(str(root / 'CellModeller/GUI/PyGLGUI.ui'))
        viewer = ui.PyGLCMViewer
        viewer.setPixelRatio(app.devicePixelRatio())
        ui.label.setTextFormat(Qt.RichText)
        ui.label.setAlignment(Qt.AlignJustify)
        ui.show()
        try:
            viewer.loadModelFile(str(root / 'Examples/multigpu_stress.py'))
        except Exception as exc:
            QMessageBox.critical(ui, 'Stress test initialization failed',
                                 str(exc) + '\n\nTry a smaller target or initial population.')
            ui.close()
            return 1
        if viewer.sim is None:  # Device selection was cancelled.
            ui.close()
            return 0

        dashboard = QDialog(ui)
        dashboard.setWindowTitle('Multi-GPU stress test — live status')
        box = QVBoxLayout(dashboard)
        status = QLabel()
        status.setTextFormat(Qt.PlainText)
        status.setTextInteractionFlags(Qt.TextSelectableByMouse)
        box.addWidget(status)
        explanation = QLabel('Work counts are cumulative; scratch is the latest stage estimate.\n'
                             'These are not GPU utilization or whole-device VRAM measurements.\n'
                             'Use Run in the colony window to start or pause.\n'
                             'Status refreshes between simulation steps; long steps can delay it.')
        box.addWidget(explanation)

        def refresh():
            sim = viewer.sim
            if sim is None:
                return
            cfg = getattr(sim.module, '_cfg', {})
            count, target = len(sim.cellStates), cfg.get('max_cells', sim.phys.max_cells)
            lines = ['Step: %d    Cells: %d / %d' % (sim.stepNum, count, target),
                     'Target reached — growth stopped' if count >= target else 'Growing toward target',
                     'Memory mode: ' + sim.clMultiGPUMemory, '']
            work = sim.CLWorkStats.get('physics.find_contacts')
            memory = sim.CLMemoryStats.get('physics.find_contacts', [])
            for index, device in enumerate(sim.CLDevices):
                lines.append('GPU %d: %s — weight %.1f%%' %
                             (sim.clDeviceNums[index], device.name, 100 * sim.clDeviceWeights[index]))
                lines.append('  Contact work: %s    Planned contact scratch: %s' %
                             (str(work[index]) if work is not None else 'not instrumented / pending',
                              '%.1f MiB' % (memory[index]['bytes'] / 2**20)
                              if index < len(memory) else 'not instrumented / pending'))
            status.setText('\n'.join(lines))

        timer = QTimer(dashboard)
        timer.timeout.connect(refresh)
        timer.start(1000)
        refresh()
        dashboard.show()
        return app.exec_()
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


if __name__ == '__main__':
    sys.exit(main())
