# -*- coding: utf-8 -*-
"""Test forms."""

from {{cookiecutter.website_name}}.public.forms import LoginForm
from {{cookiecutter.website_name}}.user.forms import RegisterForm


class TestRegisterForm:
    """Register form."""

    def test_validate_user_already_registered(self, user):
        """Enter username that is already registered.

        :param user: 

        """
        form = RegisterForm(username=user.username, email='foo@bar.com',
                            password='example', confirm='example')

        assert form.validate() is False
        assert 'Username already registered' in form.username.errors

    def test_validate_email_already_registered(self, user):
        """Enter email that is already registered.

        :param user: 

        """
        form = RegisterForm(username='unique', email=user.email,
                            password='example', confirm='example')

        assert form.validate() is False
        assert 'Email already registered' in form.email.errors

    def test_validate_success(self, db):
        """Register with success.

        :param db: 

        """
        form = RegisterForm(username='newusername', email='new@test.test',
                            password='example', confirm='example')
        assert form.validate() is True


class TestLoginForm:
    """Login form."""

    def test_validate_success(self, user):
        """Login successful.

        :param user: 

        """
        user.set_password('example')
        user.save()
        form = LoginForm(username=user.username, password='example')
        assert form.validate() is True
        assert form.user == user

    def test_validate_unknown_username(self, db):
        """Unknown username.

        :param db: 

        """
        form = LoginForm(username='unknown', password='example')
        assert form.validate() is False
        assert 'Unknown username' in form.username.errors
        assert form.user is None

    def test_validate_invalid_password(self, user):
        """Invalid password.

        :param user: 

        """
        user.set_password('example')
        user.save()
        form = LoginForm(username=user.username, password='wrongpassword')
        assert form.validate() is False
        assert 'Invalid password' in form.password.errors

    def test_validate_inactive_user(self, user):
        """Inactive user.

        :param user: 

        """
        user.active = False
        user.set_password('example')
        user.save()
        # Correct username and password, but user is not activated
        form = LoginForm(username=user.username, password='example')
        assert form.validate() is False
        assert 'User not activated' in form.username.errors
