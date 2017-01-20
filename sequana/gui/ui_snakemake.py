# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'snakemake.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Snakemake(object):
    def setupUi(self, Snakemake):
        Snakemake.setObjectName("Snakemake")
        Snakemake.resize(426, 484)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Snakemake.sizePolicy().hasHeightForWidth())
        Snakemake.setSizePolicy(sizePolicy)
        Snakemake.setAcceptDrops(False)
        Snakemake.setAutoFillBackground(False)
        Snakemake.setSizeGripEnabled(False)
        self.gridLayout = QtWidgets.QGridLayout(Snakemake)
        self.gridLayout.setObjectName("gridLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.tabs = QtWidgets.QTabWidget(Snakemake)
        self.tabs.setObjectName("tabs")
        self.tab_local = QtWidgets.QWidget()
        self.tab_local.setObjectName("tab_local")
        self.verticalLayoutWidget_2 = QtWidgets.QWidget(self.tab_local)
        self.verticalLayoutWidget_2.setGeometry(QtCore.QRect(19, 19, 141, 81))
        self.verticalLayoutWidget_2.setObjectName("verticalLayoutWidget_2")
        self.vertical_layout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_2)
        self.vertical_layout.setContentsMargins(5, 5, 5, 5)
        self.vertical_layout.setObjectName("vertical_layout")
        self.horizontal_layout = QtWidgets.QHBoxLayout()
        self.horizontal_layout.setContentsMargins(5, 5, 5, 5)
        self.horizontal_layout.setObjectName("horizontal_layout")
        self.snakemake_options_local_cores_label = QtWidgets.QLabel(self.verticalLayoutWidget_2)
        self.snakemake_options_local_cores_label.setObjectName("snakemake_options_local_cores_label")
        self.horizontal_layout.addWidget(self.snakemake_options_local_cores_label)
        self.snakemake_options_local_cores_value = QtWidgets.QSpinBox(self.verticalLayoutWidget_2)
        self.snakemake_options_local_cores_value.setMinimum(1)
        self.snakemake_options_local_cores_value.setMaximum(100000)
        self.snakemake_options_local_cores_value.setObjectName("snakemake_options_local_cores_value")
        self.horizontal_layout.addWidget(self.snakemake_options_local_cores_value)
        self.vertical_layout.addLayout(self.horizontal_layout)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.vertical_layout.addItem(spacerItem)
        self.tabs.addTab(self.tab_local, "")
        self.tab_cluster = QtWidgets.QWidget()
        self.tab_cluster.setObjectName("tab_cluster")
        self.verticalLayoutWidget_3 = QtWidgets.QWidget(self.tab_cluster)
        self.verticalLayoutWidget_3.setGeometry(QtCore.QRect(10, 20, 169, 111))
        self.verticalLayoutWidget_3.setObjectName("verticalLayoutWidget_3")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_3)
        self.verticalLayout_2.setContentsMargins(5, 5, 5, 5)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.snakemake_options_cluster_cluster_label = QtWidgets.QLabel(self.verticalLayoutWidget_3)
        self.snakemake_options_cluster_cluster_label.setObjectName("snakemake_options_cluster_cluster_label")
        self.horizontalLayout.addWidget(self.snakemake_options_cluster_cluster_label)
        self.snakemake_options_cluster_cluster_value = QtWidgets.QLineEdit(self.verticalLayoutWidget_3)
        self.snakemake_options_cluster_cluster_value.setObjectName("snakemake_options_cluster_cluster_value")
        self.horizontalLayout.addWidget(self.snakemake_options_cluster_cluster_value)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.snakemake_options_cluster_jobs_label = QtWidgets.QLabel(self.verticalLayoutWidget_3)
        self.snakemake_options_cluster_jobs_label.setObjectName("snakemake_options_cluster_jobs_label")
        self.horizontalLayout_2.addWidget(self.snakemake_options_cluster_jobs_label)
        self.snakemake_options_cluster_jobs_value = QtWidgets.QSpinBox(self.verticalLayoutWidget_3)
        self.snakemake_options_cluster_jobs_value.setSuffix("")
        self.snakemake_options_cluster_jobs_value.setMinimum(1)
        self.snakemake_options_cluster_jobs_value.setMaximum(10000)
        self.snakemake_options_cluster_jobs_value.setObjectName("snakemake_options_cluster_jobs_value")
        self.horizontalLayout_2.addWidget(self.snakemake_options_cluster_jobs_value)
        self.verticalLayout_2.addLayout(self.horizontalLayout_2)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_2.addItem(spacerItem1)
        self.tabs.addTab(self.tab_cluster, "")
        self.tab_general = QtWidgets.QWidget()
        self.tab_general.setObjectName("tab_general")
        self.formLayoutWidget_3 = QtWidgets.QWidget(self.tab_general)
        self.formLayoutWidget_3.setGeometry(QtCore.QRect(20, 10, 331, 166))
        self.formLayoutWidget_3.setObjectName("formLayoutWidget_3")
        self.layout_general = QtWidgets.QFormLayout(self.formLayoutWidget_3)
        self.layout_general.setContentsMargins(0, 0, 0, 0)
        self.layout_general.setObjectName("layout_general")
        self.snakemake_options_general_forceall_value = QtWidgets.QCheckBox(self.formLayoutWidget_3)
        self.snakemake_options_general_forceall_value.setObjectName("snakemake_options_general_forceall_value")
        self.layout_general.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.snakemake_options_general_forceall_value)
        self.snakemake_options_general_quiet_value = QtWidgets.QCheckBox(self.formLayoutWidget_3)
        self.snakemake_options_general_quiet_value.setObjectName("snakemake_options_general_quiet_value")
        self.layout_general.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.snakemake_options_general_quiet_value)
        self.snakemake_options_general_no__hooks_value = QtWidgets.QCheckBox(self.formLayoutWidget_3)
        self.snakemake_options_general_no__hooks_value.setObjectName("snakemake_options_general_no__hooks_value")
        self.layout_general.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.snakemake_options_general_no__hooks_value)
        spacerItem2 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.layout_general.setItem(7, QtWidgets.QFormLayout.FieldRole, spacerItem2)
        self.snakemake_options_general_custom = QtWidgets.QLineEdit(self.formLayoutWidget_3)
        self.snakemake_options_general_custom.setObjectName("snakemake_options_general_custom")
        self.layout_general.setWidget(6, QtWidgets.QFormLayout.FieldRole, self.snakemake_options_general_custom)
        self.label = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.label.setObjectName("label")
        self.layout_general.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.label)
        self.snakemake_options_general_verbose_value = QtWidgets.QCheckBox(self.formLayoutWidget_3)
        self.snakemake_options_general_verbose_value.setObjectName("snakemake_options_general_verbose_value")
        self.layout_general.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.snakemake_options_general_verbose_value)
        self.snakemake_options_general_summary_value = QtWidgets.QCheckBox(self.formLayoutWidget_3)
        self.snakemake_options_general_summary_value.setObjectName("snakemake_options_general_summary_value")
        self.layout_general.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.snakemake_options_general_summary_value)
        self.tabs.addTab(self.tab_general, "")
        self.verticalLayout.addWidget(self.tabs)
        self.gridLayout.addLayout(self.verticalLayout, 0, 0, 1, 1)
        self.buttonBox = QtWidgets.QDialogButtonBox(Snakemake)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout.addWidget(self.buttonBox, 1, 0, 1, 1)

        self.retranslateUi(Snakemake)
        self.tabs.setCurrentIndex(2)
        self.buttonBox.accepted.connect(Snakemake.accept)
        self.buttonBox.rejected.connect(Snakemake.reject)
        QtCore.QMetaObject.connectSlotsByName(Snakemake)
        Snakemake.setTabOrder(self.tabs, self.snakemake_options_local_cores_value)
        Snakemake.setTabOrder(self.snakemake_options_local_cores_value, self.snakemake_options_cluster_cluster_value)
        Snakemake.setTabOrder(self.snakemake_options_cluster_cluster_value, self.snakemake_options_cluster_jobs_value)

    def retranslateUi(self, Snakemake):
        _translate = QtCore.QCoreApplication.translate
        Snakemake.setWindowTitle(_translate("Snakemake", "Snakemake options"))
        self.tabs.setToolTip(_translate("Snakemake", "Snakemake parameters related to the cluster"))
        self.snakemake_options_local_cores_label.setText(_translate("Snakemake", "cores"))
        self.tabs.setTabText(self.tabs.indexOf(self.tab_local), _translate("Snakemake", "&Local"))
        self.snakemake_options_cluster_cluster_label.setToolTip(_translate("Snakemake", "<html><head/><body><p>Execute snakemake rules with the given submit command,</p><p>e.g. qsub. Snakemake compiles jobs into scripts that</p><p>are submitted to the cluster with the given command,</p><p>once all input files for a particular job are present.</p><p>The submit command can be decorated to make it aware</p><p>of certain job properties (input, output, params,</p><p>wildcards, log, threads and dependencies (see the</p><p>argument below)), e.g.: $ snakemake --cluster \'qsub</p><p>-pe threaded {threads}\'.</p><p><br/></p></body></html>"))
        self.snakemake_options_cluster_cluster_label.setText(_translate("Snakemake", "cluster"))
        self.snakemake_options_cluster_jobs_label.setText(_translate("Snakemake", "jobs"))
        self.tabs.setTabText(self.tabs.indexOf(self.tab_cluster), _translate("Snakemake", "&Cluster"))
        self.snakemake_options_general_forceall_value.setToolTip(_translate("Snakemake", "--forceall Force the execution of the selected (or the first)\n"
"rule and all rules it is dependent on regardless of\n"
"already created output.\n"
""))
        self.snakemake_options_general_forceall_value.setText(_translate("Snakemake", "forceall"))
        self.snakemake_options_general_quiet_value.setToolTip(_translate("Snakemake", "Do not output any progress or rule information"))
        self.snakemake_options_general_quiet_value.setText(_translate("Snakemake", "Quiet"))
        self.snakemake_options_general_no__hooks_value.setToolTip(_translate("Snakemake", "--no-hooks : Do not invoke onstart, onsuccess or onerror hooks   after execution.\n"
""))
        self.snakemake_options_general_no__hooks_value.setText(_translate("Snakemake", "nohooks"))
        self.label.setToolTip(_translate("Snakemake", "<html><head/><body><p>Add any other valid snakemake options here</p></body></html>"))
        self.label.setText(_translate("Snakemake", "any other options"))
        self.snakemake_options_general_verbose_value.setToolTip(_translate("Snakemake", "<html><head/><body><p> Print debugging output.</p><p><br/></p></body></html>"))
        self.snakemake_options_general_verbose_value.setText(_translate("Snakemake", "verbose"))
        self.snakemake_options_general_summary_value.setToolTip(_translate("Snakemake", "<html><head/><body><p>Print a summary of all files created by the workflow.</p><p>The has the following columns: filename, modification time, rule version, status, plan. Thereby rule version contains the versionthe file was created with (see the version keyword of rules), and status denotes whether the file is missing, its input files are newer or if version or implementation of the rule changed since file creation. Finally the last column denotes whether the file will be updated or created during the next workflow execution.</p><p><br/></p></body></html>"))
        self.snakemake_options_general_summary_value.setText(_translate("Snakemake", "summary"))
        self.tabs.setTabText(self.tabs.indexOf(self.tab_general), _translate("Snakemake", "&General"))

