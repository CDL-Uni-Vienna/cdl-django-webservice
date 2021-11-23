from django.urls import path


from knox import views as knox_views

from cdl_rest_api import views


urlpatterns = [
    # as_view is the standard function to convert APIView class
    # to be rendered by url
    path("login/", views.UserLoginApiView.as_view()),
    path("logout/", knox_views.LogoutView.as_view()),
    path("experiments/", views.ExperimentListView.as_view()),
    path("experiments/<slug:experiment_id>/", views.ExperimentDetailView.as_view()),
    path("register/", views.RegisterView.as_view()),
    path("update/<int:pk>/", views.UpdateView.as_view()),
]
