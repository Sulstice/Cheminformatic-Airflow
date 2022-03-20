from flask import Blueprint
from flask_appbuilder import expose, BaseView as AppBuilderBaseView

from airflow.plugins_manager import AirflowPlugin
from flask_admin.base import MenuLink

bp = Blueprint(
    "plotting_plugins", __name__,
    template_folder='templates', # registers airflow/plugins/templates as a Jinja template folder
    static_folder='static',
    static_url_path='/static/plotting_plugins ')

# Creating a flask appbuilder BaseView
class PCAAnalysisAppBuilderBaseView(AppBuilderBaseView):

    template_folder = '/home/sulstice/airflow/plugins/plotting_plugins/templates'
    
    @expose("/")
    def list(self):
    
        return self.render_template("pca_analysis.html")

class ProbabilityDistributionAppBuilderBaseView(AppBuilderBaseView):

    template_folder = '/home/sulstice/airflow/plugins/plotting_plugins/templates'

    @expose("/")
    def list(self):

        return self.render_template("probability_distribution.html")

class CGenFFCompassAppBuilderBaseView(AppBuilderBaseView):

    template_folder = '/home/sulstice/airflow/plugins/plotting_plugins/templates'
    
    @expose("/")
    def list(self):

        return self.render_template("cgenff_compass.html")

class ZMatrixStoreAppBuilderBaseView(AppBuilderBaseView):

    template_folder = '/home/sulstice/airflow/plugins/plotting_plugins/templates'

    @expose("/")
    def list(self):
        
        return self.render_template("zmatrix_store.html")

pca_analysis_appbuilder_view = PCAAnalysisAppBuilderBaseView()
pca_analysis_appbuilder_package = {
        "name": "PCA Analysis",
        "category": "Pipeline Plots",
        "view": pca_analysis_appbuilder_view
}

probability_distribution_appbuilder_view = ProbabilityDistributionAppBuilderBaseView()
probability_distribution_appbuilder_package = {
        "name": "Probability Distribution Penalty Scores",
        "category": "Pipeline Plots",
        "view": probability_distribution_appbuilder_view
}

cgenff_compass_appbuilder_view = CGenFFCompassAppBuilderBaseView()
cgenff_compass_appbuilder_package = {
        "name": "CGenFF Compass",
        "category": "Pipeline Plots",
        "view": cgenff_compass_appbuilder_view
}

z_matrix_appbuilder_view = ZMatrixStoreAppBuilderBaseView()
z_matrix_appbuilder_package = {
        "name": "Z-Matrix Store",
        "category": "Pipeline Plots",
        "view": z_matrix_appbuilder_view,
}


# Defining the plugin class
class AirflowTestPlugin(AirflowPlugin):
    name = "plotting_plugins"
    # operators = []
    # flask_blueprints = [bp]
    # hooks = []
    appbuilder_views = [pca_analysis_appbuilder_package, probability_distribution_appbuilder_package, cgenff_compass_appbuilder_package, z_matrix_appbuilder_package]
    # executors = []
    # admin_views = []
