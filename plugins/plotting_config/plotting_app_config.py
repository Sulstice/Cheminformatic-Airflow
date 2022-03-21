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

pca_analysis_appbuilder_view = PCAAnalysisAppBuilderBaseView()
pca_analysis_appbuilder_package = {
        "name": "PCA Analysis",
        "category": "Pipeline Plots",
        "view": pca_analysis_appbuilder_view
}



# Defining the plugin class
class AirflowTestPlugin(AirflowPlugin):
    name = "plotting_plugins"
    # operators = []
    # flask_blueprints = [bp]
    # hooks = []
    appbuilder_views = [pca_analysis_appbuilder_package]
    # executors = []
    # admin_views = []
